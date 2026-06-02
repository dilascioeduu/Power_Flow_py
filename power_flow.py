import numpy as np
import pandas as pd
import time
import scipy.linalg as la

import os
os.system('cls' if os.name == 'nt' else 'clear')


def power_flow(arquivo_excel):

    arquivo_excel = caminho

    print('=== IMPORTANDO DADOS DO EXCEL ===')

    # Configuração das bases globais
    Sbase = float(input('Informe a potência base do sistema(MVA): \n'))
    print('=== PROCESSANDO DADOS BRUTOS NO PADRÃO IEEE ===')

    # 1. Leitura direta das abas do Excel usando pandas (equivalente ao readtable)
    tabela_raw_barras = pd.read_excel(arquivo_excel, sheet_name='Barras')
    tabela_raw_linhas = pd.read_excel(arquivo_excel, sheet_name='Linhas')

    nb = len(tabela_raw_barras)
    nl = len(tabela_raw_linhas)

    # 2. Pré-alocação da matriz final: [ID, Tipo, V, Theta, P_pu, Q_pu]
    barras = np.zeros((nb, 6))

    for i in range(nb):
        barras[i, 0] = tabela_raw_barras['Bus_No'].iloc[i]  # Copia o ID da Barra
        
        # Tradução automática do código de operação IEEE (1 - Slack; 2 - PV; 3 - PQ)
        ieee_code = tabela_raw_barras['Bus_Code'].iloc[i]
        if ieee_code == 1:
            barras[i, 1] = 1  # 1 vira Slack 
        elif ieee_code == 2:
            barras[i, 1] = 2  # 2 vira PV 
        else:
            barras[i, 1] = 3  # 0 vira PQ 
        
        barras[i, 2] = tabela_raw_barras['V_pu'].iloc[i]  # Módulo de tensão inicial
        barras[i, 3] = tabela_raw_barras['V_th'].iloc[i]  # Ângulo inicial (0 rad)
        
        # Adaptando as potências líquidas da tabela: (Geração - Carga) / Sbase
        barras[i, 4] = (tabela_raw_barras['Gen_MW'].iloc[i] - tabela_raw_barras['Load_MW'].iloc[i]) / Sbase
        barras[i, 5] = (tabela_raw_barras['Gen_Mvar'].iloc[i] - tabela_raw_barras['Load_Mvar'].iloc[i]) / Sbase

    # 3. Montagem direta da matriz de linhas 
    linhas = np.column_stack((
        tabela_raw_linhas['From_Bus'],
        tabela_raw_linhas['To_Bus'],
        tabela_raw_linhas['R_pu'],
        tabela_raw_linhas['X_pu'],
        tabela_raw_linhas['B_pu'],
        tabela_raw_linhas['Tap_pu']
    ))

    print('>> Dados convertidos para pu e condicionados com sucesso!\n')

    # Configurações de execução
    tol = float(input('\nDigite a tolerância (ex: 1e-4): '))
    max_iter = 20
    tic = time.time()  # Início do cronômetro (equivalente ao tic)

    # 3 - Formação da matriz de admitância (Y BUS)
    Ybus = np.zeros((nb, nb), dtype=complex)

    for k in range(nl):
        # np.where é o equivalente do comando find() do MATLAB, retornando os índices onde a condição é verdadeira.
        # Pega o índice base 0 da barra correspondente
        i = np.where(barras[:, 0] == linhas[k, 0])[0][0] # Resgata as barras de origem (De) que o cabo k interliga
        j = np.where(barras[:, 0] == linhas[k, 1])[0][0] # Resgata as barras de destino (Para) que o cabo k interliga
        
        z = linhas[k, 2] + 1j * linhas[k, 3]  # Monta a impedância (R + jX)
        B_shunt = linhas[k, 4]
        tap = linhas[k, 5]
        
        y = 1 / z  # Transforma em admitância
        
        Ybus[i, j] = -y / tap # A admitância é dividida pelo tap para refletir a transformação de tensão (valor negativo pois estamos fora da diagonal principal)
        Ybus[j, i] = -y / tap
        Ybus[i, i] = Ybus[i, i] + (y / (tap**2)) + 1j * (B_shunt / 2) # A admitância é dividida pelo tap ao quadrado para...
        Ybus[j, j] = Ybus[j, j] + y + 1j * (B_shunt / 2) # ...refletir a transformação de tensão, e o shunt é adicionado
        # A admitância é adicionada diretamente para a barra de destino, e o shunt é adicionado

    Ym = np.abs(Ybus) # Armazena o valor do módulo de cada elemento da matriz
    Yth = np.angle(Ybus) # Armazena o valor do ângulo de cada elemento da matriz (em radianos)

    # 4 - Início do processo iterativo (NR)
    V = np.copy(barras[:, 2]) # Inicializa o vetor de tensões com os valores iniciais da tabela
    theta = np.copy(barras[:, 3]) # Inicializa o vetor de ângulos com os valores iniciais da tabela (0 rad para todas as barras, exceto a slack que pode ter um valor diferente)
    iter_count = 0
    conv = False

    while not conv and iter_count < max_iter:
        iter_count += 1
        Pcalc = np.zeros(nb)
        Qcalc = np.zeros(nb)

        for i in range(nb):
            for j in range(nb):
                # alpha = gamaij + thetai - thetaj
                ang = Yth[i, j] + theta[j] - theta[i] # O ângulo é a soma do ângulo da admitância e a diferença dos ângulos de tensão
                Pcalc[i] += Ym[i, j] * V[i] * V[j] * np.cos(ang) # Potência ativa injetada na barra
                Qcalc[i] -= Ym[i, j] * V[i] * V[j] * np.sin(ang) # Potência reativa injetada na barra 

        dP = barras[:, 4] - Pcalc # Mismatch de potência ativa: potência líquida (Pesp) da barra (geração - carga) menos a potência calculada (Pcalc)
        dQ = barras[:, 5] - Qcalc # Mismatch de potência reativa: potência líquida (Qesp) da barra (geração - carga) menos a potência calculada (Qcalc)

        idx_P = np.where(barras[:, 1] != 1)[0] # Procura na matriz onde a barra não é referência (PQ ou PV precisam de dP)
        idx_Q = np.where(barras[:, 1] == 3)[0] # Procura na matriz onde a barra é do tipo PQ (apenas PQ precisam de dQ) 
        
        # Concatenação do vetor Mismatches deltaP e deltaQ
        mis = np.concatenate((dP[idx_P], dQ[idx_Q]))

        if np.max(np.abs(mis)) < tol: # Verifica o critério de parada
            conv = True
            break

        # Montagem da Matriz Jacobiana
        H = np.zeros((nb, nb))
        N = np.zeros((nb, nb))
        M = np.zeros((nb, nb))
        L = np.zeros((nb, nb))
        
        for i in range(nb):
            for j in range(nb):
                ang = Yth[i, j] + theta[j] - theta[i]
                if i != j: # Verifica se estamos fora da diagonal principal
                    H[i, j] = -V[i] * V[j] * Ym[i, j] * np.sin(ang)
                    N[i, j] = V[i] * Ym[i, j] * np.cos(ang)
                    M[i, j] = -V[i] * V[j] * Ym[i, j] * np.cos(ang)
                    L[i, j] = -V[i] * Ym[i, j] * np.sin(ang)
                else: # Cálculos utilizando as simplificações
                    H[i, i] = -Qcalc[i] - (V[i]**2 * Ybus[i, i].imag)
                    N[i, i] = (Pcalc[i] / V[i]) + (V[i] * Ybus[i, i].real)
                    M[i, i] = Pcalc[i] - (V[i]**2 * Ybus[i, i].real)
                    L[i, i] = (Qcalc[i] / V[i]) - (V[i] * Ybus[i, i].imag)

        # Montagem da Jacobiana combinando os quadrantes usando np.ix_ para mimetizar o MATLAB...
        # ...eliminando as linhas e colunas correspondentes às barras de referência (Slack) e às barras PV para os cálculos de dP e dQ, respectivamente
        J11 = H[np.ix_(idx_P, idx_P)] # Quadrante H da Jacobiana, selecionando apenas as linhas e colunas correspondentes às barras que não são do tipo Slack (PQ ou PV)
        J12 = N[np.ix_(idx_P, idx_Q)] # Quadrante N da Jacobiana, selecionando apenas as linhas correspondentes às barras que não são do tipo Slack (PQ ou PV) e as colunas correspondentes às barras do tipo PQ
        J21 = M[np.ix_(idx_Q, idx_P)] # Quadrante M da Jacobiana, selecionando apenas as linhas correspondentes às barras do tipo PQ e as colunas correspondentes às barras que não são do tipo Slack (PQ ou PV)
        J22 = L[np.ix_(idx_Q, idx_Q)] # Quadrante L da Jacobiana, selecionando apenas as linhas e colunas correspondentes às barras do tipo PQ
        # A função np.ix_ é utilizada para criar um índice de seleção que permite extrair submatrizes específicas da matriz original, mimetizando o comportamento do comando de indexação do MATLAB.
        
        J_sup = np.hstack((J11, J12)) # Combina os quadrantes H e N horizontalmente para formar a parte superior da Jacobiana
        J_inf = np.hstack((J21, J22)) # Combina os quadrantes M e L horizontalmente para formar a parte inferior da Jacobiana
        J = np.vstack((J_sup, J_inf)) # Combina as partes superior e inferior verticalmente para formar a matriz Jacobiana completa, pronta para a resolução do sistema linear
        # A função np.hstack é utilizada para combinar as matrizes horizontalmente, enquanto a função np.vstack é utilizada para combinar as matrizes verticalmente, mimetizando o comportamento do comando de concatenação do MATLAB.

        # Fatoração LU e resolução de estado exatamente como o "[L, U] = lu(J); U\(L\mis)"
        # A função np.linalg.solve utiliza fatoração LU nativamente para solucionar o sistema
        dx = np.linalg.solve(J, mis)

        n_p = len(idx_P)
        theta[idx_P] = theta[idx_P] + dx[0:n_p] # Atualização do estado angular para a próxima iteração
        V[idx_Q] = V[idx_Q] + dx[n_p:] # Atualização do estado de tensão para a próxima iteração

    # 5 - Cálculo dos fluxos (P e Q)
    print('\n=PROCESSAMENTO DE FLUXOS E PERDAS=')

    P_fluxo = np.zeros((nl, 2)) # [P_de_para, P_para_de]
    Q_fluxo = np.zeros((nl, 2)) # [Q_de_para, Q_para_de]
    P_perda = np.zeros(nl)
    Q_perda = np.zeros(nl)

    for k in range(nl):
        i = np.where(barras[:, 0] == linhas[k, 0])[0][0] # i = linhas(k,0)
        #  Resgata as barras de origem (De) que o cabo k interliga 
        j = np.where(barras[:, 0] == linhas[k, 1])[0][0] # j = linhas(k,1)
        # Resgata as barras de destino (Para) que o cabo k interliga
        
        y_linha = 1 / (linhas[k, 2] + 1j * linhas[k, 3]) # Adamitância do próprio cabo
        b_metade = 1j * (linhas[k, 4] / 2) # Admitância shunt do cabo, dividida por 2 para refletir a distribuição simétrica entre as barras de origem e destino
        tap = linhas[k, 5] # Tap do transformador, utilizado para ajustar a tensão e o fluxo de potência entre as barras de origem e destino
        
        Vi = V[i] * np.exp(1j * theta[i]) # Tensão complexa na barra de origem (De), representada em forma polar (magnitude e ângulo)
        Vj = V[j] * np.exp(1j * theta[j]) # Tensão complexa na barra de destino (Para), representada em forma polar (magnitude e ângulo)
        
        # Potência complexa S = V * I_conj, onde I é a corrente e pode ser expressa em função das tensões e da admitância da linha.
        S_ij = Vi * np.conj(((Vi / (tap**2)) - (Vj / tap)) * y_linha + Vi * b_metade) 
        S_ji = Vj * np.conj((Vj - (Vi / tap)) * y_linha + Vj * b_metade)
        
        P_fluxo[k, 0] = S_ij.real
        Q_fluxo[k, 0] = S_ij.imag
        P_fluxo[k, 1] = S_ji.real
        Q_fluxo[k, 1] = S_ji.imag
        
        P_perda[k] = (S_ij + S_ji).real
        Q_perda[k] = (S_ij + S_ji).imag

    Perdas_Totais_Ativas_pu = np.sum(P_perda)
    Perdas_Totais_Reativas_pu = np.sum(Q_perda)

    # Conversão para grandezas reais usando a Potência Base (MW e Mvar)

    Perdas_Totais_Ativas_MW = Perdas_Totais_Ativas_pu * Sbase
    Perdas_Totais_Reativas_Mvar = Perdas_Totais_Reativas_pu * Sbase

    # Definição das bases para conversão
    Vbase = float(input('\nDigite a tensão base do sistema(kV): '))

    Nom_kV = np.ones(nb) * Vbase
    PU_Volt = np.copy(V)
    Volt_kV = V * Vbase
    Angle_Deg = np.rad2deg(theta)

    Gen_MW = np.zeros(nb)
    Gen_Mvar = np.zeros(nb)
    Load_MW = np.zeros(nb)
    Load_MVar = np.zeros(nb)

    From_Number = linhas[:, 0].astype(int)
    To_Number = linhas[:, 1].astype(int)

    MW_From = P_fluxo[:, 0] * Sbase
    MVar_From = Q_fluxo[:, 0] * Sbase
    MVA_From = np.sqrt(MW_From**2 + MVar_From**2)

    MW_To = P_fluxo[:, 1] * Sbase
    MVar_To = Q_fluxo[:, 1] * Sbase
    MVA_To = np.sqrt(MW_To**2 + MVar_To**2)

    MW_Loss = P_perda * Sbase
    MVar_Loss = Q_perda * Sbase

    # Criação da tabela de Linhas idêntica ao table() do MATLAB
    Tabela_Linhas_PW = pd.DataFrame({
        'From': From_Number, 'To': To_Number,
        'MW_ik': MW_From, 'Mvar_ik': MVar_From, 'MVA_ik': MVA_From,
        'MW_ki': MW_To, 'MVar_ki': MVar_To, 'MVA_ki': MVA_To,
        'MW_Loss': MW_Loss, 'Mvar_Loss': MVar_Loss
    })

    for i in range(nb):
        tipo = barras[i, 1]
        if tipo == 1:
            Gen_MW[i] = Pcalc[i] * Sbase
            Gen_Mvar[i] = Qcalc[i] * Sbase
        elif tipo == 2:
            Gen_MW[i] = barras[i, 4] * Sbase
            Gen_Mvar[i] = Qcalc[i] * Sbase
        elif tipo == 3:
            Load_MW[i] = -barras[i, 4] * Sbase
            Load_MVar[i] = -barras[i, 5] * Sbase

    # 6. APRESENTAÇÃO DOS RESULTADOS
    toc = time.time() - tic  # Equivalente ao fim do toc
    print(f'\nIterações para convergência: {iter_count}')
    print(f'Tempo de simulação: {toc:.4f} segundos')

    # 1. Puxa os dados estáticos originais
    Gen_MW_raw = tabela_raw_barras['Gen_MW'].values.copy()
    Gen_Mvar_raw = tabela_raw_barras['Gen_Mvar'].values.copy()
    Load_MW_raw = tabela_raw_barras['Load_MW'].values.copy()
    Load_Mvar_raw = tabela_raw_barras['Load_Mvar'].values.copy()

    type_list = []

    for i in range(nb):
        tipo = barras[i, 1]
        
        if tipo == 1:
            Gen_MW_raw[i] = (Pcalc[i] * Sbase) + Load_MW_raw[i]
            Gen_Mvar_raw[i] = (Qcalc[i] * Sbase) + Load_Mvar_raw[i]
            type_list.append('SLACK')
            
        elif tipo == 2:
            Gen_Mvar_raw[i] = (Qcalc[i] * Sbase) + Load_Mvar_raw[i]
            type_list.append('PV')
            
        elif tipo == 3:
            type_list.append('PQ')

    Tabela_PowerWorld = pd.DataFrame({
        'Name': tabela_raw_barras['Bus_No'],
        'Type': type_list,
        'PU_Volt': PU_Volt,
        'Volt_kV': Volt_kV,
        'Angle_Deg': Angle_Deg,
        'Gen_MW': Gen_MW_raw,
        'Gen_Mvar': Gen_Mvar_raw,
        'Load_MW': Load_MW_raw,
        'Load_Mvar': Load_Mvar_raw
    })

    # Saídas identicamente formatadas
    print('\n=====================================================================================================================')
    print('                                                     BUSES ')
    print('=====================================================================================================================')
    print(Tabela_PowerWorld.to_string(index=False))

    print('\n=====================================================================================================================')
    print('                                                   POWER FLOW ')
    print('\n=====================================================================================================================')
    print(Tabela_Linhas_PW.to_string(index=False))

    print('\n=====================================================================================================================')
    print('                         RESUMO GLOBAL DE PERDAS DA REDE                        ')
    print('=====================================================================================================================')
    print(f'Perdas Ativas Totais (P_loss): {Perdas_Totais_Ativas_pu:.4f} pu ({Perdas_Totais_Ativas_MW:.2f} MW)')
    print(f'Perdas Reativas Totais (Q_loss): {Perdas_Totais_Reativas_pu:.4f} pu ({Perdas_Totais_Reativas_Mvar:.2f} Mvar)')
    print('=====================================================================================================================')

if __name__ == "__main__":
    caminho = input("Digite o nome da planilha (ex: model.xlsx): ")
    power_flow(caminho)
    input("\nPressione ENTER para fechar...")