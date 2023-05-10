from os import sep
import numpy as np
import argparse
import sys

class lattice:
	def __init__(self, S, T, J):
		self.S = S
		self.T = T
		self.J = J
		self.spins = np.ones(S)
		self.bonds = np.array([np.ones(S),np.ones(S)])

	def random_initial_condition(self, seed=0):
		S=self.S
		#uncorrelated error distribution for all physical qubits
		for i in range(S[0]):
			for j in range(S[1]):
				if np.random.random_sample() < 1/2:
					self.spins[i][j] = -1
				else:
					self.spins[i][j] = 1

	def get_energy(self): #Calculate H
		S=self.S
		energy = 0
		for i in range(S[0]):
			for j in range(S[1]):
				energy += self.bonds[0,i,j]*self.J*self.spins[i][j]*self.spins[divmod(1+i,S[0])[1]][j]
				energy += self.bonds[1,i,j]*self.J*self.spins[i][j]*self.spins[i][divmod(1+j,S[0])[1]]
		return energy

	def print(self):
		S=self.S
		print_out = []
		# for i in range(2*S[0]):
		# 	a_line = []
		# 	for j in range(2*S[1]):
		# 		a_line.append(" ")
		# 		print_out.append(a_line)
		print(self.spins)
		print(3)

	def make_step(self):
		S = self.S
		x = np.random.randint(S[0])
		y = np.random.randint(S[1])
		d_energy = 0
		d_energy += self.bonds[0,x,y]*self.spins[x][y]*self.spins[divmod(x+1,S[0])[1]][y]
		d_energy += self.bonds[0,divmod(x-1,S[0])[1],y]*self.spins[x][y]*self.spins[divmod(x-1,S[0])[1]][y]
		d_energy += self.bonds[1,x,y]*self.spins[x][y]*self.spins[x][divmod(y+1,S[1])[1]]
		d_energy += self.bonds[1,x,divmod(y-1,S[1])[1]]*self.spins[x][y]*self.spins[x][divmod(y-1,S[1])[1]]
		d_energy *= -2*self.J
	
		if d_energy <= 0:
			self.spins[x][y] *= -1
			self.energy += d_energy
		else: 
			r = np.random.random()
			if r < np.exp(-d_energy/self.T):
				self.spins[x][y] *= -1
				self.energy += d_energy

def test_random_simulation(S, T, J, p, N, n):
	testlattice = lattice(S, T, J)
	testlattice.random_initial_condition()
	testlattice.bonds = np.round(np.random.rand(2,100,100)-p+1/2)*2-1 #p\in(0,1) probability for -1
	print(testlattice.bonds)
	testlattice.energy = testlattice.get_energy()

	energy = np.zeros(N)
	magnetization = np.zeros(N)
	for i in range(N):
		if i % (N/10) == 0:
			print(".",flush=True,end='')
		for steps in range(n):
				testlattice.make_step()	
		energy[i] = testlattice.energy
		magnetization[i] = np.sum(testlattice.spins)

	return energy, magnetization

	
def main(argv):
	parser = argparse.ArgumentParser(description='Simulation parser')

	parser.add_argument('-f','--fname',default = "test")
	parser.add_argument('--Sx',type=int,default=100)
	parser.add_argument('--Sy',type=int,default=100)
	parser.add_argument('-R','--repeat',type=int,default=1)
	parser.add_argument('-N','--steps',type=int,default=1000)
	parser.add_argument('-n','--stepsize',type=int,default=1)
	parser.add_argument('-p','--probability',type=float,default=0)
	parser.add_argument("-v", "--verbose",type=int,default=0)
	parser.add_argument("--test", action="store_true")

	parser.add_argument("--T",type=float,default=5)
	parser.add_argument("-J",type=int,default=1)
	parser.add_argument("--seed",type=int,default=0)

	parser.parse_args()
	args = parser.parse_args()

	fname = args.fname
	S = [args.Sx,args.Sy]
	R = args.repeat
	N = args.steps
	n = args.stepsize
	p = args.probability
	verbose = args.verbose
	T = args.T
	J = args.J
	seed = args.seed
	
	print("T",T,"p",p,"J",J)
	
	Ts = np.linspace(0.01,T)
	gather = np.zeros([2,len(Ts),R])
	for i in range(len(Ts)):
		for rpt in range(R):
			print("i=",i,"rpt=",rpt)
			T = Ts[i]
			#along nishimori line
			p = 1/2 * (1 + np.tanh(1/T))
			print("T",T,"p", p, "J",J)
			print("steps", N, "stepsize",n,"seed",seed)
			[energy, magnetization] = test_random_simulation(S, T, J, p, N, n)
			energy = np.flip(energy)
			magnetization = np.flip(magnetization)
			gather[0,i,rpt] = energy[0]
			gather[1,i,rpt] = magnetization[0]
	print(gather)

	np.savetxt(''.join([fname,'_e_gather.out']), gather[0].astype(int), fmt='%i')
	np.savetxt(''.join([fname,'_m_gather.out']), gather[1].astype(int), fmt='%i')

if __name__ == "__main__":
   main(sys.argv[1:])
