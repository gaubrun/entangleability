#This script verifies equation (7) in the Supplemental Material

#Given vectors x_0 x_1 x_plus x_minus in a vector space V and linear functionals f_0 f_1 f_plus f_minus in V*
#satisfying condtions (i) and (ii) from Theorem 3,
#we can define a 4x4 matrix a_{ij} = f_i(x_j) for i,j in the ordered set  {0,1,plus,minus}
#It is elementary to check that this matrix can be written as follows
#for some real numers a,b,c,s,t

#  2s+2a     0    s+a+b+t s+a-b-t
#    0     2s-2a  s-a+b-t s-a-b+t
# s+a+b+u s-a+b-u  2s+2b     0
# s+a-b-u s-a-b+u    0     2s-2b

s1=var('s1')
a1=var('a1')
b1=var('b1')
t1=var('t1')
u1=var('u1')
s2=var('s1')
a2=var('a2')
b2=var('b2')
t2=var('t2')
u2=var('u2')

#We introduce a matrix A1 corresponding to the vectors x_j and the linear forms f_i
A1=matrix([[2*s1+2*a1,0,s1+a1+b1+t1,s1+a1-b1-t1],[0,2*s1-2*a1,s1-a1+b1-t1,s1-a1-b1+t1],[s1+a1+b1+u1,s1-a1+b1-u1,2*s1+2*b1,0],[s1+a1-b1-u1,s1-a1-b1+u1,0,2*s1-2*b1]])

#We introduce a matrix A2 corresponding to the vectors y_l and the linear forms g_k
A2=matrix([[2*s2+2*a2,0,s2+a2+b2+t2,s2+a2-b2-t2],[0,2*s2-2*a2,s2-a2+b2-t2,s2-a2-b2+t2],[s2+a2+b2+u2,s2-a2+b2-u2,2*s2+2*b2,0],[s2+a2-b2-u2,s2-a2-b2+u2,0,2*s2-2*b2]])

#We introduce a vector omega corresponding to the linear combination of "x_j tensor y_l" defined in (4)
omega=vector([0,0,1,0,0,1,0,0,1,0,-1,0,0,0,0,0])

#We introduce a vector phi corresponding to the linear combination of "f_i tensor g_k" defined in (6)
phi=vector([1,3,-1,1,3,1,1,-1,-1,1,1,-1,1,-1,-1,1])

#We compute phi(omega) as a linear combination of the coefficients of the matrix "A1 tensor A2"
factor(phi*A1.tensor_product(A2)*omega)

#The output is 4*(t1 - u1)*(t2 - u2), which coincides with
#4(f_0(x_plus)-f_plus(x_0))(g_0(y_plus)-g_plus(y_0))
