代码说明

* 所有的二维代码中V必须是方阵
* 通过 solveXXX 和 eigXXX 求得的向量中表示的是数值解在基函数下的系数 需要用 getvalXX 转换成函数值
* 如何使用代码和画图可以参考 test1d 和 test2d 和 test_pnas
* test_pnas 在Neumann边界下重复了PNAS文章中的前几个图


function [U, lam] = eigD1d(V, num, N)
求解一维Dirichlet边界特征值问题

输入:
V(一维向量):   分片常数的数值
num(整数):    所求的特征值个数
N(整数):      分片多项式的次数(默认为10，如非必要，无需手动设置)

输出:
U(二维向量):   大小为(?, num) 每一列代表一个特征函数
lam(一维向量): 大小为(num, 1) 所求的特征值


function [U, lam] = eigD2d(V, num, N)
求解二维Dirichlet边界特征值问题 输入和输出类似


function [U, lam] = eigR1d(V, h, num, N)
求解一维Robin边界特征值问题

输入:
h(常数):   方程中的参数h

其它输入和输出类似


function [U, lam] = eigR2d(V, h, num, N)
求解二维Robin边界特征值问题 输入和输出类似


function U = solveD1d(V, N)
求解一维Dirichlet边界右端项问题 右端项为1 边界为0

输入:
V(一维向量):   分片常数的数值
N(整数):      分片多项式的次数(默认为10，如非必要，无需手动设置)

输出:
U(一维向量):   求出的数值解 需要用getval1d转换成函数值


function U = solveD2d(V, N)
求解二维Dirichlet边界右端项问题 输入和输出类似


function U = solveR1d(V, g, x0, u0, N)
求解一维Neumann边界右端项问题

输入:
g(常数):      方程中的参数g 外法向导数
x0, u0(常数):  强制边界条件 u(x0)=u0 默认不输入 也就是没有内部的强制 必须满足(0<x0<1)

其它输入和输出类似


function U = solveR2d(V, g, x0, y0, u0, N)
求解二维Neumann边界右端项问题 输入和输出类似


function u = getval1d(U, Ns, N)
把数值解在基函数上的投影转化成函数值

输入:
U(一维向量): 数值解在基函数上的投影
Ns(常数):    每个区间内采样点的个数 默认50
N(整数):     分片多项式的次数 必须和求解方程时的次数匹配 默认10

输出:
u(一维向量): 数值解在采样点上的函数值


function u = getval2d(U, Ns, N)
把数值解在基函数上的投影转化成函数值

输入:
U(一维向量): 数值解在基函数上的投影
Ns(常数):    每个区间内采样点的个数 默认20

输出:
u(二维向量): 数值解在采样点上的函数值

其它输入和输出类似


function u = my_nmlz(u)
归一化 把u归一化到绝对值最大是1


function [vx, vy, vw] = valley_line(w)
求出valleyline上的点

输入:
w(二维向量): 数值解在采样点上的函数值

输出:
vx,vy(一维向量): valleyline 上面所有点和横坐标和纵坐标
vw(一维向量): valleyline 上面所有点对应的函数值 画 effective valleyline 的时候会用到