
import time

def main_function(n):

    print(n)

    for i in range(n):
        print(i**n)

    print("[INFO] - Finished at time {}".format(time.time()))
