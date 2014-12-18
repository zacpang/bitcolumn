import random
f = open('/Users/zhifei/Desktop/test.txt','w')
start = 1;
end = pow(2,32)-1;
for i in range(1,1000000):
    a =int(random.uniform(start, end));
    b =int(random.uniform(start, end));
    c =int(random.uniform(start, end));
    d =int(random.uniform(start, end));
    e =int(random.uniform(start, end));
    f.write(str(a)+","+str(b)+","+str(c)+","+str(d)+","+str(e))
    f.write('\n')
f.close();
