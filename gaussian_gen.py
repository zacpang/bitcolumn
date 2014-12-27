import random
out = open("10attr_data.txt","w");
for i in range(1,25600000):
        a = 0;b = 0;c = 0;d = 0;e = 0; f = 0; g = 0; h = 0;j = 0; k = 0;
        while(a <= 0):
                a =int(random.gauss(pow(2,20),pow(2,25)));
        while(b <= 0):
                b =int(random.gauss(pow(2,20),pow(2,25)));
        while(c <= 0):
                c =int(random.gauss(pow(2,20),pow(2,25)));
        while(d <= 0):
                d =int(random.gauss(pow(2,20),pow(2,25)));
        while(e <= 0):
                e =int(random.gauss(pow(2,20),pow(2,25)));
        while(f <= 0):
                f =int(random.gauss(pow(2,20),pow(2,25)));
        while(g <= 0):
                g =int(random.gauss(pow(2,20),pow(2,25)));
        while(h <= 0):
                h =int(random.gauss(pow(2,20),pow(2,25)));
        while(j <= 0):
                j =int(random.gauss(pow(2,20),pow(2,25)));
        while(k <= 0):
                k =int(random.gauss(pow(2,20),pow(2,25)));
        out.write(str(a)+','+str(b)+','+str(c)+','+str(d)+','+str(e)+','+str(f)+','+str(g)+','+str(h)+','+str(j)+','+str(k));
        out.write('\n');
out.close();
