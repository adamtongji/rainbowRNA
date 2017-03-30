import os,sys
from lib.decorator import time_dec,time_func
import time


@time_dec
def main():
    time.sleep(1)



main()

