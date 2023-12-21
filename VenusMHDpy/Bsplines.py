import numpy as np

def alphaBspline(x,i,j,grid,d=0):

	R = np.zeros_like(x)
	a = grid.knots[i+grid.sk]
	b = grid.knots[i+j+grid.sk]
	c = abs(a - b) >= 1.E-12
	if d == 0: R[c] = (x[c]-a[c])/(b[c]-a[c])
	if d == 1: R[c] = 1./(b[c]-a[c])

	return R

def Bspline(x,i,j,grid,der=0):
	f = np.zeros_like(x)

	if j==0:
		if der == 0: f[(grid.knots[i+grid.sk] <= x) & (x < grid.knots[i+grid.sk+1]) | (i == grid.N+j) & (x == grid.S[-1])] = 1.
		if der == 1: f = np.zeros_like(x)
	else:
		if der == 0: f = alphaBspline(x,i,j,grid)*Bspline(x,i,j-1,grid)+(1.-alphaBspline(x,i+1,j,grid))*Bspline(x,i+1,j-1,grid)
		if der == 1: f = j*(alphaBspline(x,i,j,grid,d=1)*Bspline(x,i,j-1,grid)-alphaBspline(x,i+1,j,grid,d=1)*Bspline(x,i+1,j-1,grid))

	# Bspline j=2. I implemented this to check if the recursion Bspline function above was slow.
	#It turns out that this implementation is faster, but just a little bit. Not fast enough so that I 
	#sacrifice simplicity of the code vs performance
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#B_0i = np.zeros_like(x)
	#B_0ip1 = np.zeros_like(x)
	#B_0ip2 = np.zeros_like(x)

	#B_0i[(grid.knots[i+grid.sk] <= x) & (x < grid.knots[i+grid.sk+1]) | (i == grid.N+j) & (x == grid.S[-1])] = 1.
	#B_0ip1[(grid.knots[i+1+grid.sk] <= x) & (x < grid.knots[i+1+grid.sk+1]) | (i+1 == grid.N+j) & (x == grid.S[-1])] = 1.
	#B_0ip2[(grid.knots[i+2+grid.sk] <= x) & (x < grid.knots[i+2+grid.sk+1]) | (i+2 == grid.N+j) & (x == grid.S[-1])] = 1.
	
	#alphar = alphaBspline(x,i+1,1,grid)

	#f = alphaBspline(x,i,2,grid)*(alphaBspline(x,i,1,grid)*B_0i+(1.-alphar*B_0ip1))+(1.-alphaBspline(x,i+1,2,grid))*(alphar*B_0ip1+(1.-alphaBspline(x,i+2,1,grid)*B_0ip2))
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	return f





