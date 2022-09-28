# -*- coding: utf-8 -*-
"""
AUTOVALORES E AUTOVETORES DE MATRIZES TRIDIAGONAIS SIMÉTRICAS - O ALGORITMO QR

- Virgínia Boina Dalvi.

"""

import numpy as np            #Chamada para a biblioteca NumPy (trabalho com vetores e matrizes)

#-----------------------------------FUNÇÕES-----------------------------------#

# Função que retorna a matriz nxn cuja diagonal principal vale 2 e as imediatas subdiagonais abaixo e acima valem -1:

def A(n):
   A = np.zeros([n,n])
   np.fill_diagonal(A,2)
   np.fill_diagonal(A[:-1, 1:],-1)
   np.fill_diagonal(A[1:, :-1],-1)
   return (A)

# Função que retorna o seno e cosseno a ser utilizado na matriz de Givens:

def sc(matriz, i):
    if abs(matriz[i][i])>abs(matriz[i+1][i]):
        tau = matriz[i+1][i]/matriz[i][i]
        c = 1/np.sqrt(1+(np.square(tau)))
        s = c*tau
        return(s,c)
    else:
        tau = matriz[i][i]/matriz[i+1][i]
        s = 1/np.sqrt(1+(np.square(tau)))
        c = s*tau
        return(s,c)

# Função que monta a matriz de Givens:
    
def givens(n,i,j,s,c):
    G = np.identity(n)
    G[i][i] = c
    G[j][j] = c
    G[i][j] = s
    G[j][i] = -s
    return(G)

# Função que executa a rotação de Givens, retorna a matriz resultante triangular superior R e guarda os 
# valores de seno e cosseno de cada iteração:

def R(matriz,n):
    list = []
    for k in range(n-1):
        ang = sc(matriz,k)
        list.append(ang)
        G = givens(n,k,k+1,ang[0],ang[1])
        matriz = np.matmul(G,matriz)
    return(matriz,list)

# Função que executa o algoritmo QR e devolve a matriz de autovalores A e a matriz 
# de seus respectivos autovetores V, com a precisão desejada:
    
def QR(matriz,n,e):
    simetrica = matriz
    autovalores = R(matriz,n)[0]
    autovetores = np.identity(n)
    count = 0
    m = 0
    while m>=0:
        for k in range(n-1):
            ang = R(simetrica,n)[1][k]
            autovalores = np.matmul(autovalores,np.transpose(givens(n,k,k+1,ang[0],ang[1])))
            autovetores = np.matmul(autovetores,np.transpose(givens(n,k,k+1,ang[0],ang[1])))
        count = count + 1
        guarda = autovalores
        simetrica = guarda
        for all in range(n-1):
            if abs(simetrica[k+1][k])<e:
                return(autovalores,autovetores,count)
            else:
                autovalores = R(guarda,n)[0]

# Função que retorna os autovalores e seus autovetores associados com a precisão desejada:
    
def auto(matriz,n,e):
    listaA = []
    listaV = []
    count2 = 0
    for i in range(n,1,-1):
        matrizA = QR(matriz,i,e)[0]
        matrizV = QR(matriz,i,e)[1]
        count2 = count2 + 1
        listaA.append(matrizA[i-1][i-1])
        listaV.append(matrizV[i-1][i-1])
        matriz = matrizA[:i-1,:i-1]
        if i == 2:
            listaA.append(matrizA[i-2][i-2])
            listaV.append(matrizV[i-2][i-2])
    return(listaA,listaV,count2)

#----------------------------------- MAIN ------------------------------------#

def main():
    
    print("----------------")
    print("BEM-VINDO!!!!")
    print("----------------")
    print("Aqui vão as soluções!")
    print("")
    print("Iremos aplicar o método da fatoração QR para as matrizes tridiagonais simétricas com a diagonal principal igual a 2 e as subdiagonais iguais a -1:")
    print()
    print("Para n = 2, temos:")
    print()
    print("Matriz inicial:")
    print()
    print(A(2))
    print()
    print("Autovalores e autovetores calculados pelo método QR com a precisão de 0.000001:")
    print()
    print("Os autovalores são:")
    print()
    print(auto(A(2),2,0.000001)[0])
    print()
    print("Os autovetores associados são:")
    print()
    print(auto(A(2),2,0.0000001)[1])
    print()
    print("O número de iterações do médoto foi:")
    print(QR(A(2),2,0.000001)[2])
    print()
    print("Para n = 4, temos:")
    print()
    print("Matriz inicial:")
    print()
    print(A(4))
    print()
    print("Autovalores e autovetores calculados pelo método QR com a precisão de 0.000001:")
    print()
    print("Os autovalores são:")
    print()
    print(auto(A(4),4,0.000001)[0])
    print()
    print("Os autovetores associados são:")
    print()
    print(auto(A(4),4,0.0000001)[1])
    print()
    print("O número de iterações do médoto foi:")
    print(QR(A(4),4,0.000001)[2])
    print()
    print("Para n = 8, temos:")
    print()
    print("Matriz inicial:")
    print()
    print(A(8))
    print()
    print("Autovalores e autovetores calculados pelo método QR com a precisão de 0.000001:")
    print()
    print("Os autovalores são:")
    print()
    print(auto(A(8),8,0.000001)[0])
    print()
    print("Os autovetores associados são:")
    print()
    print(auto(A(8),8,0.0000001)[1])
    print()
    print("O número de iterações do médoto foi:")
    print(QR(A(8),8,0.000001)[2])
    print()
    print("Para n = 16, temos:")
    print()
    print("Matriz inicial:")
    print()
    print(A(16))
    print()
    print("Autovalores e autovetores calculados pelo método QR com a precisão de 0.000001:")
    print()
    print("Os autovalores são:")
    print()
    print(auto(A(16),16,0.000001)[0])
    print()
    print("Os autovetores associados são:")
    print()
    print(auto(A(16),16,0.0000001)[1])
    print()
    print("O número de iterações do médoto foi:")
    print(QR(A(16),16,0.000001)[2])
    print()
    print("Para n = 32, temos:")
    print()
    print("Matriz inicial:")
    print()
    print(A(32))
    print()
    print("Autovalores e autovetores calculados pelo método QR com a precisão de 0.000001:")
    print()
    print("Os autovalores são:")
    print()
    print(auto(A(32),32,0.000001)[0])
    print()
    print("Os autovetores associados são:")
    print()
    print(auto(A(32),32,0.0000001)[1])
    print()
    print("O número de iterações do médoto foi:")
    print(QR(A(32),32,0.000001)[2])
    print()
    print()
    print('Fim!')

main()