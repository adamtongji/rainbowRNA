#!/usr/bin/env python
#coding:utf-8
import os,sys,time
import collections
from functools import wraps


def const():
    """
    this function defines inner constants, which will not change
    in whole scripts
        Usage:
                import const (class _const saved as const.py)
                const.ConstNames = ConstValues
        Notice:
                ConstNames should be all upperCase!
                this func(class) can save as another file for import

    """
    class _const:
        class ConstError(TypeError):
            pass

        class ConstCaseError(ConstError):
            pass

        def __setattr__(self, name, value):
            if self.__dict__.has_key(name):
                raise self.ConstError, "Can't change const.%s" % name
            if not name.isupper():
                raise self.ConstCaseError,\
                    'const name "%s" is not all uppercase!' % name
            self.__dict__[name] = value

    sys.modules[__name__] = _const()


def time_dec(func):
    @wraps(func)
    def _time_dec(*args, **kwargs):
        const.TITLE = 'RNAseq pipeline'
        const.AUTHOR = os.popen("whoami").read().rstrip()
        const.TEMPLATEDATE = 'Dec 22, 2016'
        const.CODEFUNCTION = 'This program is for rnaseq.'
        begin_time = time.clock()
        print "\033[1;31;38m"
        print "#" * 50
        print 'User ' + const.AUTHOR + ':'
        print const.TITLE + ' runs at ' + \
              time.strftime("%Y-%m-%d %X", time.localtime()) + '...'
        print 'Pipeline update: ' + const.TEMPLATEDATE
        print const.CODEFUNCTION
        print "#" * 50
        print "\033[0m"
        res = func(*args, **kwargs)
        end_time = time.clock()
        lapsed_time = end_time - begin_time

        print "\033[1;31;38m"
        print "#" * 50
        print const.TITLE + ' ends at ' + \
              time.strftime("%Y-%m-%d %X", time.localtime()) + '...'
        print "Totally \t%.03f seconds lapsed" % lapsed_time
        print "#" * 50
        print "\033[0m"

        return res

    return _time_dec


def time_func(func):
    @wraps(func)
    def _time_func(*args, **kwargs):
        start= time.clock()
        res = func(*args, **kwargs)
        end = time.clock()
        lapse = end - start
        print "{} takes {} seconds.".format(func.__name__, lapse)
        return res
    return _time_func


# incompleted
def stat_log(logfile):
    def _stat_log(func):
        @wraps(func)
        def __stat_log(*args, **kwargs):
            with open(logfile,"a") as f:
                pass