class Chess(object):
    def __init__(self, k):
        self.k = k
        self.k_p = 2**k
        self.board = [([0] * self.k_p) for i in range(self.k_p)]
        self.xpos = -1 + 2**(self.k-1)
        self.ypos = -1 + 2**(self.k-1)
        self.xbou = 0
        self.ybou = 0

    def output(self):
        for i in range(self.k_p):
            for j in range(self.k_p):
                print(self.board[i][j], end='\t\t')
            print('\n')

class L(Chess):
    count = 0

    def filling(self, k, x, y, xdir, ydir):
        while(k > 1):
            self.filling(k-1, x - ((-1)**(1-xdir))*(2**(k-2)) + (-1)**(xdir), y + ((-1)**(1-ydir))*(2**(k-2)), 1-xdir, ydir)
            self.filling(k-1, x + ((-1)**(1-xdir))*(2**(k-2)), y - ((-1)**(1-ydir))*(2**(k-2)) + (-1)**(ydir), xdir, 1-ydir)
            self.filling(k-1, x + ((-1)**(1-xdir))*(2**(k-2)), y + ((-1)**(1-ydir))*(2**(k-2)), xdir, ydir)
            k -= 1
        self.board[x][y] = self.count + self.k*(10**(self.k))
        self.board[x + (-1)**xdir][y] = self.count + self.k*(10**(self.k))
        self.board[x][y + (-1)**ydir] = self.count + self.k*(10**(self.k))
        self.count += 1

class P(Chess):
    def setting(self, x, y):
        self.x = x
        self.y = y
        self.board[x][y] = 0

    def adding_1(self):
        l = L(self.k)
        l.filling(self.k, int(2**(self.k-1)), int(2**(self.k-1)), 1, 1)
        for i in range(l.k_p):
            for j in range(l.k_p):
                xb=((self.xpos+(self.x <= self.xpos)) > self.x)
                yb=((self.ypos+(self.y <= self.ypos)) > self.y)
                self.board[self.xbou+i][self.ybou+j] += l.board[ (-1)**(1-xb)*i + (1-xb)*(l.k_p-1) ][ (-1)**(1-yb)*j + (1-yb)*(l.k_p-1) ]
        self.xpos += ((-1)**(self.x <= self.xpos)) * (2**(self.k-2))
        self.ypos += ((-1)**(self.y <= self.ypos)) * (2**(self.k-2))
        self.xbou = self.xbou + (2**(self.k-1)) * (self.x >= (self.xbou + 2**(self.k-1)))
        self.ybou = self.ybou + (2**(self.k-1)) * (self.y >= (self.ybou + 2**(self.k-1)))

    def adding(self):
        while(self.k > 0):
            self.adding_1()
            self.k -= 1



print('请输入k：')
k = input()
b=P(int(k))

print("请输入x：")
x = input()

print("请输入y：")
y = input()

b.setting(int(x), int(y))
b.adding()
b.output()