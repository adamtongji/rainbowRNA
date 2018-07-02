#!/usr/bin/env python
#coding:utf-8
import os,sys


class ConfigParser(object):

    def __init__(self):
        self.__configs = {}
        self.__preset = {"mapping":{"Treat":None,"Control":None,"HISAT2_path":None,
                                  "HISAT2_index":None,"Seqtype":None,"Inputcheck":None
                                  ,"Outputdir":None,"Max_process":None,"Annotationfile":None,
                                  "feature_count_path":None,"Genomefile":None,"ReadLength":None,
                                    "STARpath":None, "STARindex":None, "STARindexdir":None},
                       "deseq":{"Pair_rep":None,"pvalue":None,"Outputdir":None,
                                "Expr_dir":None},
                       "downstream":{"Outputdir":None,"Genome":None},
                       "circ_mapping":{"Treat":None,"Control":None,"Seqtype":None,
                                       "Inputcheck":None,"Outputdir":None,"Genomefile":None,
                                       "Max_process":None,"Annotationfile":None},
                       "circ_process":{"Outputdir":None,"Pair_rep":None,"pvalue":None,
                                          "Genomefile": None,"Genome":None,"Treat":None,
                                          "Annotationfile":None}
        }

    def __runtype_config(self,_runtype):
        if len(_runtype)<1:
            print "Unexpected error in runtype!"
            sys.exit(1)
        self.runtype = _runtype
        for _type in self.runtype:
            self.__configs[_type] = {}

    def __single_config(self, _runtype,param, val):
        if not _runtype in self.runtype:
            print "Omit {} paramters: {}".format(_runtype,param)
        else:
            self.__configs[_runtype][param] = val
            print "{} paramters {} is {}".format(_runtype,param,val)

    def __multi_config(self, _runtype, param, val):
        if not _runtype in self.runtype:
            print "Omit {} paramters: {}".format(_runtype, param)
        else:
            self.__configs[_runtype][param] = val.split(",")
            print "{} paramters {} is {}".format(_runtype, param, val)

    # dump is used for debugging
    def load_configs(self,configfile="config.txt"):
        _file = [i.strip() for i in open(configfile)]
        for _line in _file:
            if ":" in _line[:40] and _line.split(":")[0]=="Runtype":
                mode = _line.split(":")[1].lower()
                if mode == "full":
                    runtype = ("mapping", "deseq", "downstream")
                elif mode =="de_full":
                    runtype = ("deseq", "downstream")
                elif mode == "mapping":
                    runtype = ("mapping")
                elif mode == "de_only":
                    runtype = ("deseq")
                elif mode == "circ_full":
                    runtype = ("circ_mapping","circ_process")
                elif "circ" in mode:
                    runtype = tuple([mode])
                else:
                    print "Error in runtype!"
                    sys.exit(1)
        self.__runtype_config(runtype)
        for _line in _file:
            if ":" in _line[:40]:
                _key = _line.split(":")[0]
                _value = _line.split(":")[1]
                for _type in runtype:
                    if self.__preset[_type].has_key(_key):
                        if "," in _value:
                            self.__multi_config(_type,_key,_value)
                        else:
                            self.__single_config(_type,_key,_value)
        Configs = self.__configs
        for _type in runtype:
            for _key,_val in Configs[_type].iteritems():
                if _val is None:
                    Configs[_type].pop(_key)

        return runtype, Configs


def input_check(*files):
    for _files in files[0]:
        if _files.endswith("gz"):
            # print "If your run ciri,please uncompress the gzip file ^_^"
            pass
            # sys.exit(1)
        elif _files.endswith("fastq") or _files.endswith("fq"):
            continue
        else:
            print _files
            print "Unknown input format"
            sys.exit(1)

    _file0 = files[0]+files[1]
    if len(_file0)>len(set(_file0)):
            print "Same files in input.  Please check treat and control filenames!"
            sys.exit(1)

    for _file in files:
        for each in _file:
            if not os.path.isfile(each):
                print '{0} do not exist!'.format(each)
                sys.exit(1)
