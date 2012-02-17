#!/usr/bin/env python

import os

# unix only workaround for making sure to capture ALL stdout to the gui

if __name__ == "__main__":
    os.system("python interface.pyw 2>&1 | tee output.log")
