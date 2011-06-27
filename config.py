from pyplasm import *
from os import getcwd


COMPOUNDS_DIR = "".join([getcwd(), "/compounds/"])


if __name__ == '__main__':
    print COMPOUNDS_DIR
    print type(COMPOUNDS_DIR)