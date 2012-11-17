import os
import sys
import sqlite3

from subprocess import call
from pkg_resources import resource_filename

#from chilin.MotifParser import MotifParser

from jinja2 import Environment, PackageLoader

exists = os.path.exists

class PipeController(object):
    def __init__(self, conf, rule, log, **args):
        """
        read in Options from command line
        Get template and conf information
        """
        self.cmd = ''
        self.shellextract = 0
        self.env = Environment(loader=PackageLoader('chilin', 'template'))
        self.conf = conf
        self.rule = rule
        self.log = log
        self.debug = args.get("debug", False)
        self.threads = args.get("threads", 1)
        self.datasummary = args.get("datasummary", 1)

    def run_cmd(self, cmd, exit_ = True, error_handler = lambda :False):
        """
        univeral call shell and judge
        """
        self.log("Run command:\t"+cmd)
        if call(cmd, shell = True):
            # if encounters error
            result = error_handler()
            if exit_:
                print "`"+cmd+"`"+" failed, exit"
                sys.exit(0)
            else:
                return result
        else:
            return True

    def ifnot_runcmd(self, test, cmd, exit_ = True, else_handler = lambda :True, size_check=True):
        """
        if the test is string and the file path of it exists, skip this cmd
        else if the test is NOT string and it's True, skip this cmd
        else execute the cmd
        """
        not_null = lambda x: os.path.isfile(x) and os.path.getsize(x) >0


        if type(test) == str:
            self.log("checking file " + test)
            if not exists(test):
                test = False     
            else:
                if not_null(test) and size_check:
                    test = True
                else:
                    test = False
        if not test:
            return self.run_cmd(cmd, exit_)
        else:
            else_handler()
            return self.log(cmd+" is skipped")

    def smart_run(self, cmd, can_skip = False, can_exit = True):
        if self.debug:
            self.ifnot_runcmd(can_skip, cmd, can_exit)
        else: self.run_cmd(cmd, can_exit)

    def cp(self, orig, new):
        return 'cp -rf %s %s' % (orig, new)

    def _render(self):
        """
        write into the DA.txt template
        """
        DA = self.env.get_template('DA.txt')
        ds_rendered = DA.render(self.rendercontent)
        with open(self.datasummary,"a") as df:
            df.write(ds_rendered)
            
class QC_Controller(object):
    def __init__(self, conf, rule, log, texfile, **args):
        """
        All the class in the module derives from this class
        """

        jinja_env = Environment(loader = PackageLoader('chilin', 'template'),
                                block_start_string = '\BLOCK{',
                                block_end_string = '}',
                                variable_start_string = '\VAR{',
                                variable_end_string = '}',
                                comment_start_string = '\#{',
                                comment_end_string = '}',
                                line_statement_prefix = '%-',
                                line_comment_prefix = '%#',
                                trim_blocks = True,
                                autoescape = False,
                                )
        self.env = jinja_env
        self.render = {}
        self.checker = []
        self.summaryCheck = []
        self.summaryRender = {}
        self.conf = conf
        self.rule = rule
        self.texfile = texfile
        self.log = log
        self.debug = args.get("debug", False)
        self.threads = args.get("threads", 1)
        self.template = self.env.get_template('template.tex')
        print os.listdir(resource_filename("chilin", "template"))
        print 1
        self.Rtemplate = self.env.get_template('Rtemplate.tex')
        self.db = sqlite3.connect(resource_filename('chilin', 'db/QC_based_infomation.db')).cursor()
        self.shiftsize = args.get("shiftsize", "")

    def run_cmd(self, cmd, exit_ = True, error_handler = lambda :True):
        """
        universal call shell and judge
        """
        self.log("Run command:\t"+cmd)
        if call(cmd, shell = True):
            # If running command encounters error,
            # call() returns a non-zero value

            result = error_handler()
            if exit_:
                raise
                sys.exit(0)
            else:
                return result
        else:
            return True

    def if_runcmd(self, condition, cmd, else_handler = lambda :True, size_check = True):
        """
        run a command conditionally
        """

        if type(condition) == str:
            self.log("checking file "+condition)
            if not exists(condition):
                condition = True
            else:
                if os.path.isfile(condition) and os.path.getsize(condition) <=0 and size_check:
                    condition = True
                else:
                    condition = False

        if condition:
            return self.run_cmd(cmd)
        else:
            else_handler()
            return self.log(cmd+" is skipped")
        
    
    def _check(self):
        """ Check whether the quality of the dataset is ok.
        vcb is the callback function for value"""
        ord = lambda x:[x["desc"], x["data"], x["value"], x["cutoff"], x["test"]]
        if len(self.checker)!=0:
            for i in self.checker:
                value = round(float(i["value"]),3)
                cutoff = round(float(i["cutoff"]),3)
                if value >= cutoff:
                    i["test"] = 'Pass'
                    self.summarycheck.append(ord(i))
                else:
                    i["test"] = 'Fail'
                    self.summarycheck.append(ord(i))

    def _render(self, mode="a"):
        """ Generate the latex code for current section. """
        content = self.template.render(self.render).replace('%','\\%')
        with open(self.texfile, mode) as f:
            f.write(content)

    def _end(self, test = lambda : True):
        def add():
            return '\\end{document}'
        return add()

