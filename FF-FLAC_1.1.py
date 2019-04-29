# -*- coding: latin-1 -*-
import numpy  as np
import timeit
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D

#LOG
#versao 1.01 - Adicionada interpolacao nearest para pontos fora do grid de interpolacao
#versao 1.02 - Adicionado funcao teste e removidos os nan da tabela - bugfix
#versao 1.03 - Saida somente id-pp
#versao 1.1 - Implementado achamatch2D-transiente(implementado rotina arruma_nan), 
#             Arrumado achamatch3D-transiente(implementado rotina arruma_nan), 
#             Arrumados graficos 2D nao transiente ->Colocada colorbar,
#             Arrumados graficos 2D transiente ->Colocada colorbar,
#             Arrumados graficos 3D transiente com colorbar & tamanho certo, 
#             Arrumada saida id-pp para todos -> deletado outras formas de entrada
#             Conferencia geral dos resultados

def le_feeflow(input_fee_flow): #leitura da grid para lista de listas a partir de arquivo
    
    print ('Lendo linhas do flac...')
    
    from_file = open(input_fee_flow, 'r')
    linhas = from_file.readlines()
    linhas_split = []
    for linha in linhas:
        linhas_split.append(linha.split())
    header_parser = linhas_split[0]
    valores_input = linhas_split[1:]
    print ('Concluido...')
    from_file.close()
    return valores_input, header_parser
    
def le_flacgrid(input_flac_grid_coord):   #leitura da grid do FLAC para lista de listas a partir de arquivo
    print ('Lendo linhas do feeflow...')
    from_file = open(input_flac_grid_coord, 'r')
    linhas = from_file.readlines()
    linhas_split = []
    for linha in linhas:
        linha = linha.replace("(", "")
        linha = linha.replace(",", " ")
        linha = linha.replace(")", "")
        linha = linha.split()
        linhas_split.append(linha)
    
    coord_points = linhas_split
    print ('Concluido...')
    return coord_points

def acha_match_2D(entrada_ff, entrada_grid_flac, transient, address_list):#rotina de interpolacao: 2D
    if transient == False:  #1/2 - rotina de interpolacao: 2D - Nao Transiente
        
        # Identificacao das colunas do FF na lista de listas
        x_ff_address = address_list[0] 
        y_ff_address = address_list[1]
        z_ff_address = address_list[2]
        value_address = address_list[3]
        print ('Inicio da rotina de calculo de pontos.Isto pode levar algum tempo...')
        #Transformacao em numpy array
        array_ff = np.array(entrada_ff, dtype = float)
        
        array_flac = np.array(entrada_grid_flac, dtype = float)
        
        #Input de tolerancia para consideracao de pontos FF como pontos FLAC no eixo 
        tolerance = float(raw_input("""Selecione uma tolerancia para generalizacao,  
na mesma unidade que o modelo!
Ans:"""))
        #Conversao dos valores de z para mesma referencia do FLAC
        zs_flac = np.unique(array_flac[:,3])*-1
        rows_flac = array_flac.shape[0]
        id_table = np.zeros(shape=(rows_flac,4))
    
        for z in zs_flac:
            
            array_ff_filt = array_ff[(array_ff[:,z_ff_address] < z+tolerance) & 
                                     (array_ff[:,z_ff_address] > z-tolerance)] 
            
            x_ff_filt = array_ff_filt[:,x_ff_address]
            y_ff_filt = array_ff_filt[:,y_ff_address]
            values  = array_ff_filt[:,value_address]
            
            array_flac_filt = array_flac[array_flac[:,3] == z*-1] 
            
            x_flac_filt = array_flac_filt[:,1]
            y_flac_filt = array_flac_filt[:,2]
            z_flac_filt = array_flac_filt[:,3]*-1
            
            grid_general = griddata((x_ff_filt, y_ff_filt), values, (x_flac_filt, y_flac_filt), method='linear').T
            codes = array_flac_filt[:,0]
            
            if z == zs_flac[0]:
                id_table = np.array([codes,x_flac_filt, y_flac_filt,z_flac_filt,grid_general]).T
            else:
                new_table = np.array([codes,x_flac_filt, y_flac_filt,z_flac_filt,grid_general]).T
                id_table = np.vstack((id_table,new_table))
        prev_col_number = id_table.shape[0]
        
        flaw_table = (id_table[np.isnan(id_table).any(axis=1)])
        id_table = (id_table[~np.isnan(id_table).any(axis=1)])
                
        x_flaw = flaw_table[:,1]
        y_flaw =  flaw_table[:,2]
        z_flaw =  flaw_table[:,3]
        
        flaw_grid = griddata((x_ff_filt, y_ff_filt), values, (x_flaw, y_flaw), method='nearest').T
        codes = flaw_table[:,0]
        new_table = np.array([codes,x_flaw, y_flaw,z_flaw,flaw_grid]).T
        
        id_table = np.vstack((id_table,flaw_table))
        id_table = id_table[~np.isnan(id_table).any(axis=1)]
        
        new_col_number = id_table.shape[0]
        loss = prev_col_number - new_col_number
        original_ff = array_ff.shape[0]
        return id_table, loss, original_ff 

    elif transient == True: #####CALCULO 2D - TRANSIENTE#####
        total_loss = 0
        array_ff = np.array(entrada_ff, dtype = float)
        array_flac = np.array(entrada_grid_flac, dtype = float)
        
        x_ff_address = address_list[0]
        y_ff_address = address_list[1]
        z_ff_address = address_list[2]
        value_address = address_list[3]
        time_address = address_list[4]
        
        x_flac = array_flac[:,1]
        y_flac = array_flac[:,2]
        
        #Input de tolerancia para consideracao de pontos FF como pontos FLAC no eixo 
        tolerance = float(raw_input("""Selecione uma tolerancia para generalizacao,  
na mesma unidade que o modelo!
Ans:"""))
        
        #Conversao dos valores de z para mesma referencia do FLAC
        zs_flac = np.unique(array_flac[:,3])*-1
        rows_flac = array_flac.shape[0]
        id_table = np.zeros(shape=(rows_flac,4))
        
        time = array_ff[:,time_address]
        timesteps = np.unique(time)
        
        for timestep in timesteps.astype(float):
            print ('Calculando valor para o tempo %.1f' %(timestep))
            
            for z in zs_flac:
            
                array_ff_filt_z = array_ff[(array_ff[:,z_ff_address] < z+tolerance) & 
                                     (array_ff[:,z_ff_address] > z-tolerance)] 
                
                array_ff_filt = array_ff_filt_z[(array_ff_filt_z[:,time_address] == timestep)]
            
                x_ff_filt = array_ff_filt[:,x_ff_address]
                y_ff_filt = array_ff_filt[:,y_ff_address]
                values  = array_ff_filt[:,value_address]
                
                ###VERIFICACO DA ENTRADA DO FF#####
                #fig, ax = plt.subplots()  
                #ax.scatter(x_ff_filt,y_ff_filt, c=values, cmap='coolwarm')
                #ax.set_title("Entrada do feeflow no tempo %.2f - z = %.2f" % (timestep, z))
                #plt.show()
            
                array_flac_filt = array_flac[array_flac[:,3] == z*-1] 
                
                x_flac_filt = array_flac_filt[:,1]
                y_flac_filt = array_flac_filt[:,2]
                z_flac_filt = array_flac_filt[:,3]*-1
                nflac = array_flac_filt.shape[0]
                timestep_list = np.zeros(nflac) + timestep
                
                grid_general = griddata((x_ff_filt, y_ff_filt), values, (x_flac_filt, y_flac_filt), method='linear').T
                codes = array_flac_filt[:,0]
            
                if z == zs_flac[0] and timestep == timesteps[0]:
                    id_table = np.array([codes,grid_general, timestep_list, x_flac_filt, y_flac_filt,z_flac_filt]).T
                else:
                    new_table = np.array([codes,grid_general,timestep_list, x_flac_filt, y_flac_filt,z_flac_filt]).T
                    id_table = np.vstack((id_table,new_table))
                
                prev_col_number = id_table.shape[0]
        
                flaw_table = (id_table[np.isnan(id_table).any(axis=1)]) #filtra valores com falhas
                id_table = (id_table[~np.isnan(id_table).any(axis=1)]) # tira nan da id_table
                
                x_flaw = flaw_table[:,3]
                y_flaw =  flaw_table[:,4]
                z_flaw =  flaw_table[:,5]
                codes = flaw_table[:,0]
                nflaw = flaw_table.shape[0]
                timestep_list = np.zeros(nflaw) + timestep
                
                flaw_grid = griddata((x_ff_filt, y_ff_filt), values, (x_flaw, y_flaw), method='nearest').T
                new_table = np.array([codes,flaw_grid, timestep_list,x_flaw, y_flaw,z_flaw,]).T
        
                id_table = np.vstack((id_table,new_table))
        
                new_col_number = id_table.shape[0]
                
                loss = prev_col_number - new_col_number
                total_loss += total_loss + loss
                original_ff = array_ff.shape[0]
                
        return id_table, loss, original_ff 
            

def acha_match_3D(entrada_ff, entrada_grid_flac, transient, address_list):#OK
    print ('Inicio da rotina de calculo de pontos.Isto pode levar algum tempo...')
    array_ff = np.array(entrada_ff, dtype = float)
    array_flac = np.array(entrada_grid_flac, dtype = float)
    if transient == True:  ###############################rotina de interpolacao: 3D - Transiente
        x_ff_address = address_list[0]
        y_ff_address = address_list[1]
        z_ff_address = address_list[2]
        value_address = address_list[3]
        time_address = address_list[4]
        
        x_flac = array_flac[:,1]
        y_flac = array_flac[:,2]
        z_flac = np.abs(array_flac[:,3])
        nflac = array_flac.shape[0]
        
        time = array_ff[:,time_address]
        timesteps = np.unique(time)
        
        for timestep in timesteps:
            print ('Calculando valor para o passo de tempo %.2f' %(timestep))
            
            filtered = array_ff[array_ff[:,7] == timestep]
            x_ff = filtered[:,x_ff_address]
            y_ff = filtered[:,y_ff_address]
            z_ff = filtered[:,z_ff_address]
            timestep_list = np.zeros(nflac) + timestep
            
            values  = filtered[:,value_address]        
        
            grid_int = griddata((x_ff, y_ff, z_ff), values, (x_flac, y_flac, z_flac), method='linear')
            codes = array_flac[:,0]
                        
            if timestep == timesteps[0]:
                id_table = np.array([codes,grid_int, timestep_list, x_flac,y_flac,z_flac]).T
            else:
                new_table = np.array([codes,grid_int, timestep_list, x_flac,y_flac,z_flac]).T
                id_table = np.vstack((id_table, new_table))
            
            prev_col_number = id_table.shape[0]
        
            if len(np.unique(np.isnan(id_table))) != 1: # aparece true se tiver nan na lista 
                flaw_table = (id_table[np.isnan(id_table).any(axis=1)])
                id_table = id_table[~np.isnan(id_table).any(axis=1)]
                
                x_flaw = flaw_table[:,3]
                y_flaw =  flaw_table[:,4]
                z_flaw =  flaw_table[:,5]
        
                flaw_grid = griddata((x_ff, y_ff, z_ff), values, (x_flaw, y_flaw, z_flaw), method='nearest').T
                codes = flaw_table[:,0]
                nflaw = flaw_table.shape[0]
                timestep_list = np.zeros(nflaw) + timestep
                
                new_table = np.array([codes,flaw_grid,timestep_list,x_flaw, y_flaw,z_flaw]).T
                
                id_table = np.vstack((id_table,new_table))
                           
                new_col_number = id_table.shape[0]
                loss = prev_col_number - new_col_number
                original_ff = array_ff.shape[0]
        return id_table,loss, original_ff  
    
    elif transient == False:     ################################rotina de interpolacao: 3D - Nao transiente
        x_ff_address = address_list[0]
        y_ff_address = address_list[1]
        z_ff_address = address_list[2]
        value_address = address_list[3]
        
        x_ff = array_ff[:,x_ff_address]
        y_ff = array_ff[:,y_ff_address]
        z_ff = array_ff[:,z_ff_address] 
        
        values  = array_ff[:,value_address]
        
        x_flac = array_flac[:,1]
        y_flac = array_flac[:,2]
        z_flac = np.abs(array_flac[:,3]) # inverting z axis
        grid_int = griddata((x_ff, y_ff, z_ff), values, (x_flac, y_flac, z_flac), method='linear')
        codes = array_flac[:,0]
        id_table = np.array([codes,grid_int,x_flac,y_flac,z_flac]).T
        prev_col_number = id_table.shape[0]
        
        if len(np.unique(np.isnan(id_table))) == 2: # aparece true se tiver nan na lista 
            flaw_table = (id_table[np.isnan(id_table).any(axis=1)])
            id_table = (id_table[~np.isnan(id_table).any(axis=1)])
                
            x_flaw = flaw_table[:,2]
            y_flaw =  flaw_table[:,3]
            z_flaw =  flaw_table[:,4]
        
            flaw_grid = griddata((x_ff, y_ff, z_ff), values, (x_flaw, y_flaw, z_flaw), method='nearest').T
            codes = flaw_table[:,0]
            new_table = np.array([codes,flaw_grid,x_flaw, y_flaw,z_flaw]).T
        
            id_table = np.vstack((id_table,flaw_table))
            id_table = id_table[~np.isnan(id_table).any(axis=1)]
            
        new_col_number = id_table.shape[0]
        loss = prev_col_number - new_col_number
        original_ff = array_ff.shape[0]
        
        return id_table,loss, original_ff    
    
def escreve_saida_2D(id_table, output_file, runtime, ndim, transient, loss,original_ff):#SAIDA-2D
    if transient == True: 
        print ('Escrevendo arquivo de saida...')
        
        to_file = open(output_file, 'w')

        for line in range(id_table.shape[0]):
            #x = id_table[line,1]
            #y = id_table[line,2]
            #tolerance = 0.0001
            pcode = id_table[line,0]
            valor = id_table[line,1]
            tempo = id_table[line,2]
            
            to_file.write('%d %.5f timestep %.4f\n' % (pcode, valor, tempo))
            if line % 80 == 0:
                progressBar(line, id_table.shape[0])
    
    elif transient == False: #OK
        id_table = id_table[id_table[:,0].argsort()]
        print ('Escrevendo arquivo de saida...')
        to_file = open(output_file, 'w')

        for line in range(id_table.shape[0]):
            valor = id_table[line,4]
            
            pcode = id_table[line,0]
            to_file.write ('%d %.5f\n' % (pcode, valor)) #descongelar para saida normal
            if line % 80 == 0:
                progressBar(line, id_table.shape[0])
    
    to_file.close()
    print ('\nConcluido. A saida esta em:')
    print output_file
    print ('A conversao levou %.2f segundos ou %3f horas' % (runtime, runtime/3600))
    print ('Numero de pontos fornecidos pelo FF: %d' % (original_ff))
    print ('Numero de pontos utilizados pela interpolacao: %d' % (id_table.shape[0]))
    print ('Pontos fornecidos pelo FLAC: %.d' % (loss+id_table.shape[0]))
    print ('Pontos fora do alcance da interpolacao linear: %.d' % (loss))

    post_menu(id_table,ndim,transient)
    
def escreve_saida_3D(id_table, output_file, runtime, ndim, transient, loss, original_ff):#SAIDA-3D
    
   if transient == False:
       print ('Escrevendo arquivo de saida...')
       
       to_file = open(output_file, 'w')

       for line in range(id_table.shape[0]):
           pcode = id_table[line,0]
           valor = id_table[line,1]           
           
           to_file.write ('%d %.5f\n ' % (pcode, valor)) # congelar para saida alternativa
           
           if line % 80 == 0:
                progressBar(line, id_table.shape[0])
                
       to_file.close()
       print ('\nConcluido. A saida esta em:')
       print output_file
       print ('A conversao levou %.2f segundos ou %3f horas' % (runtime, runtime/3600))
       print ('Numero de pontos fornecidos pelo FF: %d' % (original_ff))
       print ('Numero de pontos utilizados pela interpolacao: %d' % (id_table.shape[0]))
       print ('Pontos fornecidos pelo FLAC: %.d' % (loss+id_table.shape[0]))
       print ('Pontos fora do alcance da interpolacao linear: %.d' % (loss))
       
       post_menu(id_table,ndim,transient)
   elif transient == True:
       print ('Escrevendo arquivo de saida...')
       
       to_file = open(output_file, 'w')

       for line in range(id_table.shape[0]):
           pcode = id_table[line,0]
           valor = id_table[line,1]
           time = id_table[line,2]

           to_file.write('initial pp %.5f range id %d timestep %.2f\n' % (valor,pcode, time))
           if line % 80 == 0:
                progressBar(line, id_table.shape[0])
                
       to_file.close()
       print ('\nConcluido. A saida esta em:')
       print output_file
       print ('A conversao levou %.2f segundos ou %3f horas' % (runtime, runtime/3600))
       print ('Numero de pontos fornecidos pelo FF: %d' % (original_ff))
       print ('Numero de pontos utilizados pela interpolacao: %d' % (id_table.shape[0]))
       print ('Pontos fornecidos pelo FLAC: %.d' % (loss+id_table.shape[0]))
       print ('Pontos fora do alcance da interpolacao linear: %.d' % (loss))
       
       post_menu(id_table,ndim,transient)

def main_menu():#Menu principal
    print """
    Transformador FeeFlow - FLAC 1.0-Python 
    Desenvolvimento - Geotecnia: Eng. Marina Trevizolli - marinatrevizolli@hotmail.com
                    - Python: Eng. David Bispo - davidbispo@hotmail.com
                    
    Se voce esta usando este programa, voce leu o documento "Leia-me"
    fornecido com o programa, e esta ciente dos termos e condicoes que 
    se aplicam ao mesmo.
    """
    select_ff = raw_input("""
Insira o endereco do arquivo *.dat do Feeflow.
Ans: """)                          
    select_ff = select_ff.replace("\\", "\\\\")
    
    select_flac = raw_input("""
Insira o endereco do arquivo *.dat ou texto do Flac.
Ans: """)
    select_flac = select_flac.replace("\\", "\\\\")
    select_out = raw_input("""
Insira o endereco do arquivo de saida para leitura pelo FISH.
Ans: """)      
                           
    transient = raw_input("""
A Analise e transiente? S/N
Ans: """)      
                           
    if transient == 'S':
        transient = True
    elif transient == 'N':
        transient = False

    ndim = int(raw_input("""
A Analise e de: 
    
  2  : 2 Dimensoes generalizada
  3  : 3 Dimensoes
    
  Ans: """ ))                       
                       
    return select_ff, select_flac, select_out, transient, ndim

def learquivos(ff_adress, flac_address, transient): #Parser para o header do FF - Pega os enderecos baseados em strings
    fee_flow_grid, ff_header  = le_feeflow(ff_adress)
    flac_grid = le_flacgrid(flac_address)
    x_ff_address = ff_header.index("X")
    y_ff_address = ff_header.index("Y")
    z_ff_address = ff_header.index("Z")
    ff_value = ff_header.index("PINIT")
    
    if  transient == True:
        try:
            ff_time = ff_header.index("Time")
            address_list = [x_ff_address, y_ff_address, z_ff_address, ff_value, ff_time]
        except:
            print "seu arquivo ff nao tem uma coluna tempo. encerrando o programa.."
            exit()
    else:
        address_list = [x_ff_address, y_ff_address, z_ff_address, ff_value]
    return fee_flow_grid, flac_grid, address_list

def progressBar(value, endvalue, bar_length=20): # Barra de progresso da escrita de arquivos

        percent = float(value) / endvalue
        arrow = '-' * int(round(percent * bar_length)-1) + '>'
        spaces = ' ' * (bar_length - len(arrow))

        sys.stdout.write("\rPor Cento: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()

def grafico(id_table,ndim, transient):# Plotagem de graficos de verificacao
    if ndim == 2:
        if transient == True:
            timesteps = id_table[:,2]
            first_timestep = np.unique(timesteps).min()
            last_timestep = np.unique(timesteps).max()
            timesteps_iterate = [first_timestep,last_timestep]
            
            zs_flac = np.unique(id_table[:,5])
            z_mean = int(round(len(zs_flac)/2))
            for timestep in timesteps_iterate:
                
                plt_table = id_table[(id_table[:,5] == zs_flac[z_mean]) & (id_table[:,2] == timestep)]
                fig, ax = plt.subplots(figsize = (12,6))
                
                x = plt_table[:,3]
                y = plt_table[:,4]
                c = plt_table[:,1]

                sp = ax.scatter(x,y, c=c, cmap='coolwarm')
                plt.title("Plotagem para o passo de tempo %.3f para a secao media z = %.2f\n" %(timestep,z_mean))
                plt.colorbar(sp)
                plt.grid()    
                plt.show()
    
        elif transient == False:
            x = id_table[:,1]
            y = id_table[:,2]
            c = id_table[:,4]
            
            time = id_table[:,2]
            timesteps = np.unique(time)
            ultimot = timesteps[-1]
            primeirot = timesteps[0]
            
            id_tablep = id_table[id_table[:,2] == primeirot]
            id_tableu = id_table[id_table[:,2] == ultimot ]
                                  
            xp = id_tablep[:,3]
            
            fig, ax1 = plt.subplots(figsize = (12,6))
            sp = ax1.scatter(x,y, c=c, cmap='coolwarm')
            plt.colorbar(sp)
            plt.grid()    
            plt.show()
            
    elif ndim == 3:       
        if transient == True:
            time = id_table[:,2]
            timesteps = np.unique(time)
            ultimot = timesteps[-1]
            primeirot = timesteps[0]
            
            id_tablep = id_table[id_table[:,2] == primeirot ]
            id_tableu = id_table[id_table[:,2] == ultimot ]
                                  
            xp = id_tablep[:,3]
            yp = id_tablep[:,4]
            zp = id_tablep[:,5]
            cp = id_tablep[:,1]
            
            xu = id_tableu[:,3]
            yu = id_tableu[:,4]
            zu = id_tableu[:,5]
            cu = id_tableu[:,1]

            fig = plt.figure(figsize = (10,10))
            ax1 = fig.add_subplot(211, projection='3d')
            sp1 = ax1.scatter(xp,yp,zp, c=cp)
            ax1.set_title("Plotagem do timestep %.2f" % primeirot)
            plt.colorbar(sp1)
            
            ax2 = fig.add_subplot(212, projection='3d')
            sp2 = ax2.scatter(xu,yu,zu, c=cu)
            ax2.set_title("Plotagem do timestep %.2f" % ultimot)
            plt.colorbar(sp2)
           
            plt.show()
            
        elif transient == False:
            x = id_table[:,2]
            y = id_table[:,3]
            z = id_table[:,4]
            c = id_table[:,1]
            
            fig = plt.figure(figsize = (12,6))
            ax = fig.add_subplot(111, projection='3d')
            sp = ax.scatter(x,y,z, c=c)
            plt.colorbar(sp)
            plt.show()
            
            
def ajuda():# Ajuda
    print """
    Transformador FeeFlow - FLAC 1.0-Python 

O programa converte as poro-pressoes geradas pelo FeeFlow em arquivo .dat 
para uma uma grid pre-definida no FLAC. A sequencia de execucao e:
    
    1-Menu Principal
     
    E‰ possi­vel converter analises 2D e 3D, estacionarias ou transientes. A
sai­da e um arquivo de texto contendo o texto para entrada na interface FISH
do FLAC no caso 2D generalizado, e no caso do 3D, duas colunas separadas por 
tabulacoes contendo ID do ponto e valor de poro-pressao. 
    O menu principal pede pelo endereco dos arquivos a serem convertidos. O 
programa transforma o caminho (ex.: C:/Users/User/Documents/arquivo.txt) em um
caminho compati­vel com o python. Sendo assim:
    
    ***Nao use barras invertidas("\") para o endereco
    
    ***As opcoes de estacionariedade e numero de dimensoes possuem respostas
espeficas. Qualquer resposta estranha terminarÃ¡ o programa. 

    ***Compiladores python nao lidam bem com caminhos de arquivos longos demais

    ***Cuidado com seus arquivos. Eles serao sobrescritos
    
    2 - Modelos de arquivos
    Modelos de arquivos 2D e 3D estao incluidos no pacote do Transformador Fee-Flow - Flac  
    
    ***Arquivos de formatos diferentes nao funcionam no conversor
    
    3 - Post-Menu
    Apos a plotagem dos resultados, um pos-menu e lancado. E possivel realizar
    outra analise, plotar os resultados para verificacao, abrir esta ajuda
    ou sair do programa.
    
    
Desenvolvimento-Geotecnia: Eng. Marina Trevizolli- marinatrevizolli@hotmail.com
    
- Python: Eng. David Bispo - davidbispo@hotmail.com
    
    Os desenvolvedores do programa nao se responsabilizam pelo uso, ou pelo mau
    uso do programa. A avaliacao dos resultados cabe ao usuario, assim como se sao
    seguros, confiaveis e devidamente empregaveis na pratica da engenharia.
    
    """
def post_menu(id_table, ndim, transient): #Menu post-run
    
      post_option = int(raw_input("""
O que voce gostaria de fazer agora:
                          
1 - Compatibilizar outro modelo
2 - Visualizar uma plotagem do resultado
3 - Visualizar a ajuda
4 - Sair do programa
                          
Ans: """                          
                           ))
      if post_option == 1:
          runtime_auto()
      elif post_option ==2:
          grafico(id_table,ndim, transient)
          post_menu(id_table, ndim, transient)
      elif post_option ==3:    
          ajuda()
          post_menu(id_table, ndim, transient)
      elif post_option ==4:    
          sys.exit()
          
def runtime_auto():#Comeco do programa
         
    ff_address, flac_address, output_file, transient, ndim = main_menu()
    start = timeit.default_timer()
    ff_grid, flac_grid, address_list = learquivos(ff_address, flac_address, transient)
    
    if ndim == 2:
        id_table, loss, original_ff = acha_match_2D(ff_grid,flac_grid, transient, address_list)
        stop = timeit.default_timer()
        runtime = stop - start 
        escreve_saida_2D(id_table, output_file, runtime, ndim, transient, loss,original_ff)
    
    elif ndim == 3:
        id_table,loss, original_ff = acha_match_3D(ff_grid,flac_grid, transient, address_list)
        stop = timeit.default_timer()
        runtime = stop - start
        escreve_saida_3D(id_table, output_file, runtime, ndim, transient, loss, original_ff)

def runtime_teste(ff_address, flac_address, output_file, transient, ndim):#Comeco do teste ***sem interface grafica***
    start = timeit.default_timer()
    ff_grid, flac_grid, address_list = learquivos(ff_address, flac_address, transient)   
    if ndim == 2:
        id_table, loss, original_ff = acha_match_2D(ff_grid,flac_grid, transient, address_list)
        stop = timeit.default_timer()
        runtime = stop - start 
        escreve_saida_2D(id_table, output_file, runtime, ndim, transient, loss,original_ff)
    
    elif ndim == 3:
        id_table,loss, original_ff = acha_match_3D(ff_grid,flac_grid, transient, address_list)
        stop = timeit.default_timer()
        runtime = stop - start
        escreve_saida_3D(id_table, output_file, runtime, ndim, transient, loss, original_ff)

#Entrada para o teste - Troque pelos enderecos dos seus arquivos pra testar
#ff_address_teste = r'D:\OneDrive\old\flac3d_marina\pressao_calibrada3D.dat'
#flac_address_teste = r'D:\OneDrive\old\flac3d_marina\points_v4.txt'
#output_file_teste = r'D:\OneDrive\old\flac3d_marina\fish_3d.dat'
#transient_teste = False
#ndim_teste = 2
 
#runtime_teste(ff_address_teste, flac_address_teste, output_file_teste, transient_teste, ndim_teste)

#*****Descongele aqui e congele o de cima com um "#" para rodar com a interface grafica#
runtime_auto()