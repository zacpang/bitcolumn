import random
f = open('/Users/zhifei/Desktop/test.txt','w')
for i in range(1,100000000):
    a =int(random.gauss(pow(2,31),pow(2,28)))
    f.write(str(a))
    f.write('\n')
f.close();
