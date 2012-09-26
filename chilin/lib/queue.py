#!/usr/bin/python

"""
automatically DC pipeline on servers
get the CPU, memory and disk usages
and write in to queue process
"""

import os
import datetime
import time
import sys
import sqlite3
from ConfigParser import ConfigParser
def is_server_available():
    cmd_df_usage = "df -h | grep /mnt/Storage | awk '{print int($5)}'"
    disk = float(os.popen(cmd_df_usage).read())
    print "disk usage is: %.2f%%" % disk

    cmd_CPU_usage = "top -bn 1 | grep '^ *[1-9]' | awk '{cpuusage += $9} END {print cpuusage}'"
    cpu = int(os.popen(cmd_CPU_usage).read())
    print "number of used CPU is: %d\n" % cpu

    cmd_memory_usage = "top -bn 1 | grep '^ *[1-9]'| awk '{memusage += $10} END {print memusage}'"
    memory = float(os.popen(cmd_memory_usage).read())
    
    if disk < 90 and cpu < 15 and memory < 15:
        return True
    else:
        return False

def main():
    start = datetime.datetime.now()
    print start
    run_queue()
    time.sleep(5)

def write_queue(status, id):
    print status,'test'
    conn = sqlite3.connect('queue.db')
    curs = conn.cursor()
    sql_update = 'update queue set status="%s" where datasetid=%s' % (status,id)
    curs.execute(sql_update)
    time.sleep(5)

def read_queue(tempfile):
    """temp file from server3
    convert to queue.db
    Arguments:
    - `tempfile`: 20120802.txt
    """
    f = file(tempfile, 'r')
    conn = sqlite3.connect('queue.db')
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
        starttime TEXT,
        timeconsumed  TEXT,
        status TEXT
        )'''
                     )
    except: # to avoid repeat establish queue table
        conn.rollback()
    sql_insert = "insert into queue values(?, ?, ?, ?, ?, ?, ?, ?, ?)"
    for line in f.readlines()[1:]:
        line = line.strip().split()
        line.append(datetime.datetime.now())
        line.append('')
        line.append('')
        print line
        try:
            curs.execute(sql_insert, line)
            conn.commit()
        except: # to avoid column datasetid is not unique problem
            conn.rollback()
#    curs.close()
#    conn.close()
    cf = ConfigParser()
    cf.read('example.conf')
    sql_select = 'select * from queue where status = "" order by datasetid'
    curs.execute(sql_select)
    # update conf files
    cont = curs.fetchall()[0]
    print cont
    cf.set('UserInfo', 'User', 'testuser')
    cf.set('UserInfo', 'species', cont[3])
    cf.set('UserInfo', 'factor', cont[4])
    cf.set('UserInfo', 'datatype', cont[5])
    cf.set('UserInfo', 'treatid', cont[1])
    cf.set('UserInfo', 'controlid', cont[2])
    cf.write(open('queue1.txt', 'w'))
    return cont

def run_queue():
    """input queue.db
    update status in the
    queue.db file run
    Arguments:
    - `queue`: queue.db
    """
    tempfiles = os.listdir('queuefolder')
    status = ''
    for f in tempfiles:
        conf = read_queue(os.path.join('queuefolder', f))
        for id in conf:
            if status == '':
                server_available = is_server_available()
                if server_available == True:
                    print "server is available!\n"
                    # run pipeline
            #        command = "ChiLin.py %s -t %s" % (conf)
            #        subprocess.call(command)
                else:
                    print "server is not available!\n"
                    time.sleep(5)

            call = 0
            if not call:
                status = 'done'
            else:
                status = 'error'
            write_queue(status,conf[0][0]) # id

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me, see you!\n")
        sys.exit(0)
