import os
import time
import sys
import logging
import sqlite3
import datetime
import re

from ConfigParser import ConfigParser
from optparse import OptionParser
from subprocess import call

exists = os.path.exists
error   = logging.critical
warn    = logging.warning
adata = lambda d: d.rstrip('.txt')
pabs = os.path.abspath
path = lambda a, f: os.path.join(os.path.abspath(a), f)
name = lambda n: adata(n) + '.conf'
def log(info):
    """logname : specify autodc_date.log
    info for warning or error message
    """
    print info
    logger = logging.getLogger()
    d = datetime.datetime.now().strftime("%Y_%m_%d.log")
    handler = logging.FileHandler(d)
    formatter = logging.Formatter('%(asctime)s %(levelname)s :  %(message)s ')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.NOTSET)
    logger.info(info)

def is_server_available():
    cmd_df_usage = "df -h | grep /mnt/Storage | awk '{print int($5)}'" # for whole disk
    disk = float(os.popen(cmd_df_usage).read())
    # use 10% cpu as cpu cutoff
    cmd_CPU_usage = "top -bn 1 -u qinq | grep '^ *[1-9]' | awk '{if ($9>10) {cpuusage += 1}} END {print cpuusage}'" # for single user
    cpu = int(os.popen(cmd_CPU_usage).read())
    cmd_memory_usage = "top -bn 1 -u qinq | grep '^ *[1-9]'| awk '{memusage += $10} END {print memusage}'" # single user
    memory = float(os.popen(cmd_memory_usage).read())
    if disk < 90 and cpu < 15 and memory < 15:
        return True
    else:
        return False

def read_queue(idfile, dbname, outputd, confp):
    """temp file from server3
    convert to queue.db
    example: idfile: 1011.txt
    """
    conn = sqlite3.connect(pabs(dbname))
    curs = conn.cursor()
    try:
        curs.execute('''
        create table queue (
        datasetid  TEXT  primary key,
        treatid TEXT,
        controlid TEXT,
        species TEXT,
        factor TEXT,
        datatype TEXT,
        confpath TEXT,
        outputpath TEXT,
        starttime TEXT,
        timeconsumed  TEXT,
        status TEXT
        )''')
    except: # to avoid repeat establish queue table
        conn.rollback()
    conf_c = []
    species = lambda s: "mm9" if s.lower()=='mus' else "hg19"
    def atype(t):
        if re.match('[hH]\d*', t):
            return 'Histone'
        elif re.match('[dD][nN]ase$', t):
            return 'Dnase'
        else:
            return 'TF'
    with open(idfile, 'r') as f:
        sql_insert = "insert into queue values(%s)" % ','.join(['?']*11)
        for line in f.readlines():
            line = line.strip().split('\t')
            conf_c.append(os.path.basename(adata(idfile))) # dataset pc id
            conf_c.append(line[0]) # treat id, sep by ,
            conf_c.append(line[1]) # control id
            conf_c.append(species(line[2]))  #species
            conf_c.append(line[3])  # factor
            conf_c.append(atype(line[3])) # data type
            conf_c.append(os.path.join(confp, name(os.path.basename(adata(idfile))))) # conf paths
            conf_c.append(adata(outputd))  # output path for upload data

            conf_c.append(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")) # start time
            conf_c.append('')  # time consumed, update in the write_queue
            conf_c.append('')  # status, update in the write_queue

            try:
                curs.execute(sql_insert, conf_c) # write into autodc.db
                conn.commit() # remember to commit change
            except:
                conn.rollback()
    curs.close()
    conn.close()
    return adata(outputd)

def conf(dbname, confp, idfile, outputd):
    conn = sqlite3.connect(pabs(dbname))
    curs = conn.cursor()
    sql_select = 'select * from queue where status = "" order by datasetid'
    curs.execute(sql_select)
    confs = curs.fetchall()
    print os.path.join(confp, name(os.path.basename(adata(idfile))))
    curs.close()
    conn.close()
    return confs

def write_configs(conf, outputd, confp, idfile, datad):
    cf = ConfigParser()
    print conf
    cmd = "ChiLin.py gen -s %s " % conf[3] # species, specific setting
    log(cmd)
    try:
        call(cmd, shell=True)
        cf.read('ChiLinjinja.conf')
        fqpath = lambda f, suffix: path(datad, f + suffix)
        setdata = lambda t: ','.join(map(fqpath, t.split(','), ['.fastq']*len(t.split(',')))) # TODO, add other suffix option 
        cf.set('UserInfo', 'User', 'autodc')
        cf.set('UserInfo', 'species', conf[3])
        cf.set('UserInfo', 'factor', conf[4]) # factor
        cf.set('UserInfo', 'datatype', conf[5])
        cf.set('UserInfo', 'treatpath', setdata(conf[1]))
        cf.set('UserInfo', 'controlpath', setdata(conf[2]))
        cf.set('UserInfo', 'OutputDirectory', outputd)
        cf.set('UserInfo', 'datasetid', conf[0])
        print cf.sections()
        cf.write(open(os.path.join(confp, name(os.path.basename(adata(idfile)))), 'w'))
    except TypeError:
        log("TypeError")
        sys.exit(1)

def write_queue(oldstatus, newstatus, id, dbname):
    conn = sqlite3.connect(dbname)
    curs = conn.cursor()
    if oldstatus == '':
        sql_update = 'update queue set status="%s" where datasetid=%s' % (newstatus,id)
        try:
            curs.execute(sql_update)
            conn.commit()
        except:
            conn.rollback()
            log("auto dc has done the input data, restart after preparing other data")
    sql_select = 'select * from queue'
    curs.execute(sql_select)
    confs = curs.fetchall()
    print confs
    curs.close()
    conn.close()

def run_queue(queued, datad, outd):
    def status(s):
        pass
    # cmd_run = lambda conf, type: call("ChiLin.py run -c {0} -t {1}".format(conf, type),
    #                               shell = True)
    cmd_run = lambda conf, type: os.popen("ChiLin.py run -c {0} -t {1}".format(conf, type))

    while 1:
        queuefiles = os.listdir(pabs(queued)) # all pc files
        if len(queuefiles) == 0:
            log("No input data, please fetch from pc filling table, sleep~~~")
            time.sleep(60)
            continue
        d = datetime.datetime.now()
        conflist = d.strftime("%Y_%m_%d_confs")
        confp = pabs(conflist)
        processed = confp + "_processed"
        if not all(map(os.path.exists, map(adata, queuefiles))):
            map(lambda f: call('mkdir %s' % f, shell=True), map(adata, queuefiles))
        if not os.path.exists(confp):
            call('mkdir %s' % confp, shell=True)
        if not os.path.exists(processed):
            call('mkdir %s' % processed, shell=True)
        dbname = pabs('autodc.db')
        outputd = map(lambda f: read_queue(path(queued, f), dbname, path(outd, f), confp), queuefiles) # initialize db and read in

        # mv queuefiles to a processed directory with status
        confs = map(lambda f: conf(dbname, confp, path(queued, f), outd), queuefiles)[0]  # parse the db into tuples in list
        log("all availble now conf listed below")
        log('\n'.join('\t'.join(t) for t in confs))
        
        map(lambda f, c, o: write_configs(c, o, confp, path(queued, f), datad), queuefiles, confs, outputd)
        log('write confs in %s' % confp)
        log('begin to run ')
        if is_server_available():
            log("server is available, begin to process data")
            for c in confs:
                if c[9] == '': # status change
                    if not ('fail' or 'error') in cmd_run(c[6], c[5]).read():  # TODO, error extracting method
                        write_queue(c[9], 'done', c[0], dbname) # id
                        log("%s has done" % c[0])

                    else:
                        write_queue(c[9], 'error', c[0], dbname)
                        log("%s error" % c[0])
                    call('mv %s %s' % (os.path.join(pabs(queued),c[0] + '.txt'), processed), shell=True) # may be change
                time.sleep(3)
        else:
            log("Please wait server to get enough resources")

def main():
    usage = "usage: %prog -q <queue directory> -d <data directory> -o <output directory>"
    description = "ChIP-seq Auto Pipeline"
    optparser = OptionParser(version="%prog 1.00",description=description,usage=usage)
    # TODO, add threading ,Wed 26 Sep 2012 08:46:15 AM CST
    optparser.add_option("-q",dest="queued",type="str",
                         help = "Please input pc files directory")
    optparser.add_option("-d",dest="datad", type="str",
                         help = "Please input raw data directory")
    optparser.add_option("-o",dest="outd", type="str",
                         help = "Please directory containing output")
    (options,args) = optparser.parse_args()

    if not all([options.queued, options.datad, options.outd]):
        optparser.print_help()
        sys.exit(1)
    if not all(map(os.path.exists, [options.queued, options.datad, options.outd])):
        optparser.print_help()
        sys.exit(1)
    log('AutoDC begins to run, main process')
    run_queue(options.queued, options.datad, options.outd)
    time.sleep(5)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me, see you!\n")
        sys.exit(0)
