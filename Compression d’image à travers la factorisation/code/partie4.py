import numpy as np
import matplotlib.pyplot as plt 
from numpy.linalg import norm
from scipy.linalg import svd

img = plt.imread("img/decollage.png")

def S_matrix(S,n,m):
	A=np.zeros((n, m))
	k=min(n,m)
	for i in range(k):
		A[i][i]=S[i]
	return A

def compression(A,k):
	n=np.shape(A)[0]
	m=np.shape(A)[1]
	U,S,V=svd(A)
	S=S_matrix(S,n,m)
	for i in range (k,min(n,m)):
		S[i][i]=0
	return U@S@V

def compressed_img(img, rank):
	R = img[:,:,0]
	G = img[:,:,1]
	B = img[:,:,2]
	R = compression(R,rank)
	G = compression(G,rank)
	B = compression(B,rank)
	img_comp=np.dstack((R,G,B))
	return img_comp



img_comp=compressed_img(img, 100)

plt.imshow(img_comp)
plt.axis('off')
plt.savefig('Rang 100.png',bbox_inches='tight')
plt.show()

def effi(img, max_rank):
	d=[]
	for i in range(max_rank):
		img_comp=compressed_img(img, i)
		n=norm(img - img_comp)
		d.append(n)
	x = list(range(max_rank))
	y = d
	plt.plot(x, y)
	plt.xlabel("rang de compression k")
	plt.ylabel("distance entre les images")
	plt.savefig('Efficacite.png',bbox_inches='tight')
	plt.show()
	return 0

effi(img,500)

