---
title: RS码原理及柯西优化
date: 2023-10-12 12:00:00 +0800
categories: [数理基础]
tags: [前向纠错, 线性代数, FEC]
math: true
mermaid: true
---

## Intro
Forward Error Correction（前向纠错、FEC）技术可以用于在不可靠信道中传输信息，当信道中发生丢包时通过部分原始和对应的冗余信息恢复出全部的原始信息。举例来讲，数据块 *a* 和 *b* 生成了冗余数据块 *c* ，在接收端只收到 *a* 和 *c* 时，能够还原出原始数据块 *b*，省去了请求和响应重传的时间。
Reed–Solomon Error Correction（RS码、里所码）是一种FEC算法，假设原始和冗余信息块数量分别是 *n* 和 *k*，在接收端任意收到 *n* 个原始/冗余信息块时，通过RS解码都能还原出 *n* 个原始数据块。回到上文的例子，如果使用 *a* 和 *b* 这两个原始数据块（n=2），生成 *c* 和 *d* 两个冗余数据块（k=2），使用XOR操作时无法保证任意2个数据块到达后都能恢复出 *a* 和 *b*，而RS码可以保证。

RS码的基础是有限域，它定义了数据字符的加法和乘法操作，保证两个字符在加法/乘法后还在该有限域内；其次RS码的编码和解码都是基于矩阵乘法，编码时生成矩阵和原始数据向量相乘获得原始+冗余数据向量，解码时通过确定等效的生成矩阵，进行求逆后还原出原始数据向量；柯西优化是使用柯西矩阵优化编解码效率，对于上述两项基本原理没有改变。

本文将主要基于 *An Introduction to Galois Fields and Reed-Solomon Coding*[^Galois-RS-Intro] 和 *Optimizing Cauchy Reed-Solomon Codes for Fault-Tolerant Network Storage Applications*[^Caucy-optimize] 介绍RS码基本原理及柯西优化的相关内容。

## 有限域
### 1. Feild（域）
Field（域）是加法和乘法运算被定义的元素集合，其中的元素和定义的运算满足以下性质：
1. 两个元素进行加法和乘法运算后的结果还在域内（closure闭包性质），定义在域的加法和乘法运算都满足交换律和结合律
2. 域拥有加法元和乘法元: 0+a=a, 1*a=a
3. 定义了加法逆（减法）和乘法逆（除法）
4. 加法元$\neq$乘法元
5. ab=0 <=> a=0 or b=0

实际上在后续的介绍中，域的加法和乘法运算的闭包性质是最为重要的，其他性质能够协助一些原理推导。另外域的元素应理解为“存在”，而域内的加法/乘法运算应理解为“定义”，比如 $(a+b)\%N$ 和 $a\oplus b$都应该被理解为“加法运算”。

### 2. Finite Field（有限域）
Finite Field（有限域）即元素有限的域。常见的有限域有定义在模除运算*p*（p为质数）上的 $ℤ_{p}$，其元素为{0,1,2,...,p-1}，加法运算为 $(a+b)\ mod\ p$，乘法运算为 $ ab\ mod\ p $。

以下是$ℤ_{2}$的加法表和乘法表：

![add-multi](/assets/rs-caucy/add-multi.png){:height="300" .w-50}

### 3. Galois Field（伽罗瓦域）
Galois Field（伽罗瓦域）使用符号$GF(p^m)$表示，其元素是定义在$ℤ_{p}$上的m-1次多项式。展开后为 $a_{m-1}x^{m-1} +...+ a_{1}x^1 + a_{0}x^0$，其中系数 $a_{i}$ 取自集合{$0, 1, ..., p-1$}。在编码中我们常取p=2，此时m-1次多项式表达一个长度为m位的字，当m=8时则表示一个字节。

#### a. 加法运算
伽罗瓦域上的加法为对多项式逐位进行异或：$(x^2+1)+(x+1)+(x^2+x+1)=1$，根据异或的特性可以知道，加法的逆运算减法实际上也是逐位进行异或。

以下是 $GF(2^3)$ 的加法表：

![add-gf-2-3](/assets/rs-caucy/add-gf-2-3.png){:height="300" .w-50}


#### b. 乘法运算
伽罗瓦域$GF(2^m)$的最高次限制为m-1，但乘法的直接结果可能超过这个值。下式中为 $GF(2^3)$ 的加法运算，结果的最高位已经超过m-1=2

$5\cdot6=(x^2+1)(x^2+x)=x^4+x^3+x^2+x$

因此，需要找到不可约多项式g(x)来使得乘法的结果回到域内，使得$p(x)=a(x)b(x)$在域外，但$f(x)=p(x)\ mod\ g(x)$又回到域内，这个思想有些类似 $ℤ_{p}$ 上乘法的模除 p 操作。

$x^3+x+1$ 是 $GF(2^3)$ 上的一个不可约多项式，以下为 $5\cdot6$ 的结果，使用 g(x) 将其模除回域内的过程，最终计算的结果为 $x+1$ 即 3：

![multi-mod-gf-2-3](/assets/rs-caucy/multi-mod-gf-2-3.png){:height="400" .w-50}

通过逐个计算，可以得到 $GF(2^3)$ 的乘法表：

![multi-gf-2-3](/assets/rs-caucy/multi-gf-2-3.png){:height="300" .w-50}

#### c. 乘法函数
以下是乘法的编程实现，分别模拟了竖乘和长除的过程，需要注意的是 $GF(2^m)$ 中的加减法都是逐位进行异或操作
```c++
int gf_mult(int m, int poly, int v1, int v2) {
  // poly: low order terms of g(x)
  int prod = 0; 
  int k; 
  int mask;
  /* Multiply phase */
  for (k = 0; k < m; k++)
  {
    // 通过移位和异或操作模拟竖乘的过程
    if (v1 & 1) {
      prod ^= (v2 << k);
    }
    v1 >>= 1;
    if (v1 == 0) break;
  }

  /* Reduce phase */
  // 当m=3，mask=10000，即上一步结果可能出现的最高位
  mask = 1 << (2m-2); 
  for (k = m - 2; k >= 0; k--)
  {
    // 通过移位mask模拟g(x)的长除过程
    if (prod & mask) {
      // 直接去掉当前最高位
      prod &= ~mask;
      // 长除上部写上1，做异或减法
      prod ^= (poly << k);
    }
    mask >>= 1;
  }
}
```

#### 4. 快速乘法
伽罗瓦域中的乘法可以借助对数运算实现，即 $a \cdot b=log_{2}^{-1}(log_{2}(a)+log_{2}(b))$ 成立，其中右式的“+”是自然数加法而不是多项式加法。下表是 $GF(2^3)$ 上的对数表：

![log-gf-2-3](/assets/rs-caucy/log-gf-2-3.png){:height="100" .w-50}

其中需要注意以下问题：
1. 表中的所有单元格内都是遍历整个 $GF(2^3)$ 中所有元素尝试出来的，没有更快的生成方法
2. 表中有两个不定义的项，不定义的原因在此不做展开，可以参考[^Galois-RS-Intro]
3. $log_{2}^{-1}(7)$ 未定义：7减7后为0，$log_{2}^{-1}(0)=1 => 3\cdot6=1$
4. $log_{2}(0)$ 未定义：$0\cdot n=0 => 0$，不需要使用快速乘法

快速乘法表的生成较慢（因为需要遍历），但生成后能够大幅度提高乘法的运算速度，可以作为RS码编码的一种优化方案。




## Ref
[^Galois-RS-Intro]: [*An Introduction to Galois Fields and Reed-Solomon Coding*](https://people.computing.clemson.edu/~jmarty/papers/IntroToGaloisFieldsAndRSCoding.pdf)
[^Caucy-optimize]: [*Optimizing Cauchy Reed-Solomon Codes for Fault-Tolerant Network Storage Applications*](https://ieeexplore.ieee.org/document/1659489)