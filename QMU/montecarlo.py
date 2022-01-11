import random
import matplotlib.pyplot as plt
import numpy as np

#求圆周率近似值
def pi(total):
    count = 0
    for i in range(total):
        x = random.random() #随机生成0-1
        y = random.random()
        d1 = (x**2 + y**2) ** 0.5
        if d1 <=1 :
            count+=1
        else:
            pass
    print('pi = ', 4 * count/total)




