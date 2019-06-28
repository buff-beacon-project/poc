from sympy import *
import os
import binascii
import time



# choose safe prime near 2**nbits so that exponent e has large multiplicative order
def RSASafePrimeExponent(nbits,e):

	nbytes=nbits/8
	p=int(binascii.hexlify(os.urandom(nbytes)),16)
	while p<2**(nbits-1):
		p = 2*p
	#print "target value of p=",p

	# large prime of the order of 2**(nbits/2) used to ensure all small odd exponents have 
	# large multiplicative order
	p2=int(binascii.hexlify(os.urandom(nbytes/2)),16)
	p2=nextprime(p2)
	
	a2=p/(4*p2)
	
	
	finished=False
	while not finished:
		a2=a2-1
		p=4*a2*p2+3
		if p%12==11 and gcd(e,p-1)==1:
			if isprime(p):
				p1=(p-1)/2
				if isprime((p-1)/2):
					phi=p-1
					phiphi=p1-1
					finished=True
					#print "p=",p
					#print "phi=",phi
					#print "phiphi=",phiphi
					if pow(e,phiphi/p2,phi) != 1:
							finished=True
							
							
	#print "return"
	#print "p=",p
	#print "p1=",p1
	#print "p2=",p2
	#print "phi=",phi
	#print "phiphi=",phiphi
#	print "gcd(e,phi)=",gcd(e,phi)	
#	print "pow(e,phiphi,phi)=",pow(e,phiphi,phi)
#	print "pow(e,phiphi/p2,phi)=",pow(e,phiphi/p2,phi)
#	print "pplus1 with small factors removed=",pplus1
				
	return p;


# USE RSA PSRBG to generate sequence of 512 pseudorandom bis each second
# t=time in epoch (seconds since 00:00:00 01/01/1970)
# n = p*q where p and q are secret randomly chosen 1024 bit primes
# lam = lambda(n) = (p-1)*(q-12)/gcd(p-1,q-1) is Carmichael Lambda Totient
# e is exponent coprime to lam
# p and q are reset once per minute 
# x0 is secret randomly chosen seed less than n
# xt=pow(x0,pow(e,M*t,lam),n) is the secret pseudorandom number for time 
#     where is M=2**20 so the generator jumps about one million steps ahead each second
# zt = Sum[2**k (pow(xt,pow(e,k,lam),n) mod 2),{k,0,511}] is publicly displayed
#     512-bit  pseudorandom number for time t
#  note only one bit of each RSA pseudorandom number is used in calculation of zt
# hash(zt) is shorter number derived from zt
# hash of zt for NEXT second is displayed at time t
#
# ideally x0 has provably maximum order=lambda(n) modulo n, 
#     and e has provably large order modulo lambda(n).
#  this requires more complicated determination of primes


# calculate initial parameters p, q, lambda, etc
e=3
nbits=2048
nbytes=nbits/8
nbitsprime=nbits/2
nbytesprime=nbitsprime/8

#choose 2048 bit composite n=p*q
#choose two 1024 bit primes p and q without any special properties or checks other than 
# the exponent e must be coprime to phi(n) by requiring gcd(e,p-1)=gcd(e,q-1)=1

clock0=time.clock()
print "initializing safe primes p and q. This can take a few minutes"
p=RSASafePrimeExponent(nbitsprime,e)
print "safe prime p=",p
q=RSASafePrimeExponent(nbitsprime,e)
print "safe prime q=",q
n=p*q
lam=(p-1)*(q-1)/igcd(p-1,q-1)
print "n=",hex(n)
print "e=",e
print "lam=",lam

x0=int(binascii.hexlify(os.urandom(nbytes-1)),16)
maxorder=False
p1=(p-1)/2
q1=(q-1)/2
while not maxorder:
	x0=x0+1
	maxorder=True
	if pow(x0,lam/p1,n)==1:
		maxorder=False
	if pow(x0,lam/q1,n)==1:
		maxorder=False
	if pow(x0,lam/2,n)==1:
		maxorder=False
print "x0=",x0
clock1=time.clock()
print "time to choose safe primes =",clock1-clock0,"  seconds"



M=2**20
N=512
print "setup initial pseudorandom parameters"
t=int(time.time())
t0=t
xt0=pow(x0,pow(e,M*t0,lam),n)
zt0=0L
twok=1L
x=xt0
for k in range(0,N):
	zt0=zt0+twok*(x%2)
	x=pow(x,e,n)
	twok=2*twok
ht0=zt0%4294967296

t1=t0+1
xt1=pow(x0,pow(e,M*t1,lam),n)
zt1=0L
twok=1L
x=xt1
for k in range(0,N):
	zt1=zt1+twok*(x%2)
	x=pow(x,e,n)
	twok=2*twok
ht1=zt1%4294967296


t2=t1+1
xt2=pow(x0,pow(e,M*t2,lam),n)
zt2=0
twok=1L
x=xt2
for k in range(0,N):
	zt2=zt2+twok*(x%2)
	x=pow(x,e,n)
	twok=2*twok
ht2=zt2%4294967296


print "t0=",t0, time.asctime(time.gmtime(t0))
print "xt0=",hex(xt0)
print "zt0=",hex(zt0)
print "ht0=",hex(ht0)

print " "
print "t1=",t1, time.asctime(time.gmtime(t1))
print "xt1=",hex(xt1)
print "zt1=",hex(zt1)
print "ht1=",hex(ht1)

print " "
print "t2=",t2,time.asctime(time.gmtime(t2))
print "xt2=",hex(xt2)
print "zt2=",hex(zt2)
print "ht2=",hex(ht2)
print " "

#post t0
print "t0=",t0, time.asctime(time.gmtime(t0))
print "zt0=",hex(zt0)
print "ht1=",hex(ht1)
starttime=t0

tlastpost=t0

#reset
t0=t1
xt0=xt1
zt0=zt1
ht0=ht1

t1=t2
xt1=xt2
zt1=zt2
ht1=ht2

t2=t1+1
xt2=pow(x0,pow(e,M*t2,lam),n)
zt2=0
twok=1L
x=xt2
for k in range(0,N):
	zt2=zt2+twok*(x%2)
	x=pow(x,e,n)
	twok=2*twok
ht2=zt2%4294967296


print "start broadcasting"
count=0
t=time.time()
while True:
	t=int(time.time())
	if (t > tlastpost):

		if (t0%3600==0):
			print " "
			print "parameters used during previous hour were:"
			print "p=",hex(plast)
			print "q=",hex(qlast)
			print "n=",hex(nlast)
			print "e=",e
			print "lam=",hex(lamlast)
			print "x0=",hex(x0last)

			
		#post t0
		print " "
		print "timestamp, UTC time, zt, and lowest four bytes of next zt"
		print t0, "  ", time.strftime("%H:%M:%S %m/%d/%Y UTC", time.gmtime(t0))
		print hex(zt0)
		print hex(ht1)

		tlastpost=t0
		
		t0=t1
		xt0=xt1
		zt0=zt1
		ht0=ht1

		t1=t2
		xt1=xt2
		zt1=zt2
		ht1=ht2
		

		t2=t1+1
		xt2=pow(x0,pow(e,M*t2,lam),n)
		zt2=0
		twok=1L
		x=xt2
		for k in range(0,N):
			zt2=zt2+twok*(x&1)
			x=pow(x,e,n)
			twok=twok<<1
		ht2=zt2%4294967296
		
		if (t0%3600==3598):
			print "resetting primes before top of hour"
			plast=p
			qlast=q
			nlast=n
			lamlast=lam
			x0last=x0
			clock0=time.clock()

			p=RSASafePrimeExponent(nbitsprime,e)
			print "safe prime p=",p
			q=RSASafePrimeExponent(nbitsprime,e)
			print "safe prime q=",q
			n=p*q
			lam=(p-1)*(q-1)/igcd(p-1,q-1)
			print "n=",n
			print "e=",e
			print "lam=",lam

			x0=int(binascii.hexlify(os.urandom(nbytes-1)),16)
			maxorder=False
			while not maxorder:
				x0=x0+1
				maxorder=True
				if pow(x0,lam/p,n)==1:
					maxorder=False
				if pow(x0,lam/q,n)==1:
					maxorder=False
				if pow(x0,lam/2,n)==1:
					maxorder=False
			print "x0=",x0
			clock1=time.clock()
			print "time to reset primes =",clock1-clock0



			

	




		
