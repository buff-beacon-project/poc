from sympy import *
import os
import binascii
import time

#choose nbits long strong prime p such that:
#    p=2*a1*p1+1
#    p1=2*a2*p2+1
#    such that:
#    phi(p)=2*a1*p1
#    phi(phi(p))=phi(2*a1)*a2*a2*p2
#    a1 small enough so we can calculate its prime factors and totient 
#    so we can ensure exponent e has large multiplicative order modulo phi(n)
#    and so that we can ensure that x0 has maximum multiplicative order modulo n
#    of lambda(n) = (p-1)*(q-1)/gcd(p-1,q-1)
#  return values are list [p,p1,a1]
def RSAStrongPrime(nbits,e):

	nbytes=nbits/8
	p=int(binascii.hexlify(os.urandom(nbytes)),16)
	while p<2**(nbits-1):
		p = 2*p
	#print "p=",p
	
	p2=nextprime(int(binascii.hexlify(os.urandom(nbytes/2)),16))
	#print "p2=",p2
	
	nbytesa1=nbytes/4
	while nbytesa1 > 10:
		nbytesa1 -= 1
	a1=int(binascii.hexlify(os.urandom(nbytesa1)),16)
	p1=p/(2*a1)
	a2=(p1-1)/(2*p2)
	p1=2*a2*p2+1
	while not isprime(p1):
		a2=a2+1
		p1=2*a2*p2+1
	#print "p1=",p1
	
	
	a1=(p-1)/(2*p1)
	finished=False
	while not finished:
		a1=a1-1
		if gcd(e,2*a1)==1:
			p=2*a1*p1+1
			if isprime(p):
				phi=p-1
				phiphi=totient(2*a1)*(p1-1)
				#print "gcd(e,phi)=",gcd(e,phi)
				#print "phi=",phi
				#print "phiphi=",phiphi
				#print "pow(e,phiphi,phi)=",pow(e,phiphi,phi)
				#print "pow(e,phiphi/p2,phi)=",pow(e,phiphi/p2,phi)
				if pow(e,phiphi/p2,phi)!=1:
					pplus1=p+1
					pp=1
					for i in range(0,82025):
						pp=nextprime(pp)
						while pplus1%pp==0:
							pplus1 = pplus1/pp
					#print "pplus1 after prime factors less than 2^20 removed=",pplus1
					if pplus1 > pow(2,nbits/2):
						finished=True
	phip=p-1
	phifactors=factorint(2*a1)
	phifactors[p1]=1
	#print "p=",p
	#print "phi(p)=",phi
	#print "phi(phi(p))=",phiphi	
	#print "phifactors=",phifactors
	#print "p1=",p1
	#print "p2=",p2
	#print "a1=",a1
	#print "a2=",a2
			
	return [p,p1,a1];



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
print "initializing strong primes p and q. This can take a few seconds"
plist= RSAStrongPrime(nbitsprime,e)
p=plist[0]
p1=plist[1]
a1=plist[2]
print "strong prime p=",p
print "p-1=2*a1*p1 where p1=",p1," and a1=",a1
qlist= RSAStrongPrime(nbitsprime,e)
q=qlist[0]
q1=qlist[1]
b1=qlist[2]
print "strong prime q=",q
print "q-1=2*b1*q1 where q1=",q1," and b1=",b1
n=p*q
lam=(p-1)*(q-1)/igcd(p-1,q-1)
print "n=",n
print "e=",e
print "lam=",lam

#choose x0 with maximum multiplicative order
x0=int(binascii.hexlify(os.urandom(nbytes-1)),16)
maxorder=False
while not maxorder:
	x0=x0+1
	maxorder=True
	if pow(x0,lam/p1,n)==1:
		maxorder=False
	if pow(x0,lam/q1,n)==1:
		maxorder=False
	pfactors=factorint(2*a1)
	for pf in pfactors:
		if pow(x0,lam/pf,n)==1:
			maxorder=False
	qfactors=factorint(2*b1)
	for qf in qfactors:
		if pow(x0,lam/qf,n)==1:
			maxorder=False		
print "x0=",x0
print "pow(x0,lam,n)=",pow(x0,lam,n)
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
			print "This can take a few seconds"
			plast=p
			qlast=q
			nlast=n
			lamlast=lam
			x0last=x0
			clock0=time.clock()
			print "initializing strong primes p and q. This can take a few minutes"
			plist= RSAStrongPrime(nbitsprime,e)
			p=plist[0]
			p1=plist[1]
			a1=plist[2]
			print "strong prime p=",p
			print "p-1=2*a1*p1 where p1=",p1," and a1=",a1
			qlist= RSAStrongPrime(nbitsprime,e)
			q=qlist[0]
			q1=qlist[1]
			b1=qlist[2]
			print "strong prime q=",q
			print "q-1=2*b1*q1 where q1=",q1," and b1=",b1
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
				if pow(x0,lam/p1,n)==1:
					maxorder=False
				if pow(x0,lam/q1,n)==1:
					maxorder=False
				pfactors=factorint(2*a1)
				for pf in pfactors:
					if pow(x0,lam/pf,n)==1:
						maxorder=False
				qfactors=factorint(2*b1)
				for qf in qfactors:
					if pow(x0,lam/qf,n)==1:
						maxorder=False		
			print "x0=",x0
			print "pow(x0,lam,n)=",pow(x0,lam,n)
			clock1=time.clock()
			print "time to reset primes =",clock1-clock0



			

	
	




		
