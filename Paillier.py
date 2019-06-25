import random, math
from gmpy2 import is_prime, gcd, gcdext, lcm, powmod
import gmpy2

# 生成 z bit 大素数
rs = gmpy2.random_state(gmpy2.mpz(random.randint(0, 2 ** 31)))
def prime_generate(z):
	test =  gmpy2.mpz_urandomb(rs, z - 1) + 2 ** (z - 1)
	while not gmpy2.is_prime(test):	
		test = gmpy2.mpz_urandomb(rs, z - 1) + 2 ** (z - 1)
	return test

# 非递归的拓展欧几里得算法求 a mod b 逆元 y 满足 a * x = 1 (mod b)
def inv_mod(a, b):
	return gcdext(a, b)[1]

# paillier公钥
class public_key:
	def __init__(self, N, g):
		self.N = N
		self.g = g
		self.Nsqr = self.N * self.N
	def L(self, u):
		return (u - 1) / self.N
	def E(self, m, r):
		return (powmod(self.g, m, self.Nsqr) * powmod(r, self.N, self.Nsqr)) % self.Nsqr

# paillier私钥

class pirvate_key(object):
	def __init__(self, public_key, lambda_):
		self.public_key = public_key
		self.lambda_ = lambda_
	def L(self, u):
		return (u - 1) // self.public_key.N
	def D(self, c):
		a = self.L(powmod(c, self.lambda_, self.public_key.Nsqr))
		b = inv_mod(self.L(powmod(self.public_key.g, self.lambda_, self.public_key.Nsqr)), self.public_key.N)
		return a * b % self.public_key.N

def generate_paillier_key(bit):
	p = prime_generate(bit)
	q = prime_generate(bit)
	print("p   = ", hex(p))
	print("q   = ", hex(q))
	N = p * q
	g = N + 1
	lambda_ = lcm(p-1, q-1)
	puk = public_key(N, g)
	prk = pirvate_key(puk, lambda_)
	return puk, prk


if __name__ == '__main__':
	puk, prk = generate_paillier_key(512)
	print("N   = ", hex(puk.N))
	print("g   = ", hex(puk.g))
	print("λ(N)= ", hex(prk.lambda_))

	m1 = 123
	m2 = 234
	
	c1 = puk.E(m1, 37)
	c2 = puk.E(m2, 31)
	c1 = powmod(c1, 100, puk.Nsqr)

	decryt = prk.D(c1)
	print(decryt)