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
        
        # Tradução automática do código de operação IEEE
        ieee_code = tabela_raw_barras['Bus_Code'].iloc[i]
        if ieee_code == 1:
            barras[i, 1] = 1  # 1 vira Slack 
        elif ieee_code == 2:
            barras[i, 1] = 2  # 2 vira PV 
        else:
            barras[i, 1] = 3  # 0 vira PQ 
        
        barras[i, 2] = tabela_raw_barras['V_pu'].iloc[i]  # Módulo de tensão inicial
        barras[i, 3] = tabela_raw_barras['V_th'].iloc[i]  # Ângulo inicial
        
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
        # np.where é o equivalente exato do comando find()
        # Pega o índice base 0 da barra correspondente
        i = np.where(barras[:, 0] == linhas[k, 0])[0][0]
        j = np.where(barras[:, 0] == linhas[k, 1])[0][0]
        
        z = linhas[k, 2] + 1j * linhas[k, 3]  # Monta a impedância
        B_shunt = linhas[k, 4]
        tap = linhas[k, 5]
        
        y = 1 / z  # Transforma em admitância
        
        Ybus[i, j] = -y / tap
        Ybus[j, i] = -y / tap
        Ybus[i, i] = Ybus[i, i] + (y / (tap**2)) + 1j * (B_shunt / 2)
        Ybus[j, j] = Ybus[j, j] + y + 1j * (B_shunt / 2)

    Ym = np.abs(Ybus)
    Yth = np.angle(Ybus)

    # 4 - Início do processo iterativo (NR)
    V = np.copy(barras[:, 2])
    theta = np.copy(barras[:, 3])
    iter_count = 0
    conv = False

    while not conv and iter_count < max_iter:
        iter_count += 1
        Pcalc = np.zeros(nb)
        Qcalc = np.zeros(nb)

        for i in range(nb):
            for j in range(nb):
                ang = Yth[i, j] + theta[j] - theta[i]
                Pcalc[i] += Ym[i, j] * V[i] * V[j] * np.cos(ang)
                Qcalc[i] -= Ym[i, j] * V[i] * V[j] * np.sin(ang)

        dP = barras[:, 4] - Pcalc
        dQ = barras[:, 5] - Qcalc
        
        # Procura na matriz (idx_P pega diferentes de 1; idx_Q pega iguais a 3)
        idx_P = np.where(barras[:, 1] != 1)[0]
        idx_Q = np.where(barras[:, 1] == 3)[0]
        
        # Concatenação do vetor Mismatches
        mis = np.concatenate((dP[idx_P], dQ[idx_Q]))

        if np.max(np.abs(mis)) < tol:
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
                if i != j:
                    H[i, j] = -V[i] * V[j] * Ym[i, j] * np.sin(ang)
                    N[i, j] = V[i] * Ym[i, j] * np.cos(ang)
                    M[i, j] = -V[i] * V[j] * Ym[i, j] * np.cos(ang)
                    L[i, j] = -V[i] * Ym[i, j] * np.sin(ang)
                else:
                    H[i, i] = -Qcalc[i] - (V[i]**2 * Ybus[i, i].imag)
                    N[i, i] = (Pcalc[i] / V[i]) + (V[i] * Ybus[i, i].real)
                    M[i, i] = Pcalc[i] - (V[i]**2 * Ybus[i, i].real)
                    L[i, i] = (Qcalc[i] / V[i]) - (V[i] * Ybus[i, i].imag)

        # Montagem da Jacobiana combinando os quadrantes usando np.ix_ para mimetizar o MATLAB
        J11 = H[np.ix_(idx_P, idx_P)]
        J12 = N[np.ix_(idx_P, idx_Q)]
        J21 = M[np.ix_(idx_Q, idx_P)]
        J22 = L[np.ix_(idx_Q, idx_Q)]
        
        J_sup = np.hstack((J11, J12))
        J_inf = np.hstack((J21, J22))
        J = np.vstack((J_sup, J_inf))
        
        # Fatoração LU e resolução de estado exatamente como o "[L, U] = lu(J); U\(L\mis)"
        # A função np.linalg.solve utiliza fatoração LU nativamente para solucionar o sistema
        dx = np.linalg.solve(J, mis)

        n_p = len(idx_P)
        theta[idx_P] = theta[idx_P] + dx[0:n_p]
        V[idx_Q] = V[idx_Q] + dx[n_p:]

    # 5 - Cálculo dos fluxos (P e Q)
    print('\n=PROCESSAMENTO DE FLUXOS E PERDAS=')

    P_fluxo = np.zeros((nl, 2))
    Q_fluxo = np.zeros((nl, 2))
    P_perda = np.zeros(nl)
    Q_perda = np.zeros(nl)

    for k in range(nl):
        i = np.where(barras[:, 0] == linhas[k, 0])[0][0]
        j = np.where(barras[:, 0] == linhas[k, 1])[0][0]
        
        y_linha = 1 / (linhas[k, 2] + 1j * linhas[k, 3])
        b_metade = 1j * (linhas[k, 4] / 2)
        tap = linhas[k, 5]
        
        Vi = V[i] * np.exp(1j * theta[i])
        Vj = V[j] * np.exp(1j * theta[j])
        
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