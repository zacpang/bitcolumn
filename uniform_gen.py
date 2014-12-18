import random
f = open('/Users/zhifei/Desktop/test.txt','w')
start = 1;
end = pow(2,32)-1;
for i in range(1,100000000):
    a =int(random.uniform(start, end));
    f.write(str(a))
    f.write('\n')
f.close();
