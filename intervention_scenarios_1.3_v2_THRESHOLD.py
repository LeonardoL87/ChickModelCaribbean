#!/usr/bin/env python
# SEIR without demography (no vital dynamics)
# by Carlos J. Dommar (carlos_dot_dommar_at_ gmail_dot_com)
#                     (carlos_dot_dommar_at_ isglobal_dot_org)
########################################################################
# INTERVENTION SCRIPTS                                                 #
# INTERVENTION SCRIPTS                                                 #
# INTERVENTION SCRIPTS                                                 #
# modifications developed by Leonardo Lopez 
# (leonardorafael.lopez@isglobal.org)  to impose an artificial 
# threshold  on the ODEs solution are incorporated in this scripts 
########################################################################
# This is the SEIR epidemic with no births nor deaths.                 #
# coupled to s network model tath represents mobility of human hosts   #
# within a geographical domain (the Caribbean)                         #
########################################################################
########################################################################
# Modified by carlos j dommar - carlos.dommar@gmail.com                #
# to adapt it to a simple 2-patch model                                #
# - modification log                                                   #
#                                                                      #
# - pending udate info for version 0.6                                 #
#                                                                      #
# - v. 0.4.7:                                                          #
#       * I check that densities are used instread of absolut          #
#       population sizes                                               #
#       * Implemente a Block Random Network Model to model the         #
#       political structure in the network                             #
#       * Once parameterized the BRNM then simulate epidemics on reali #
#       -zations drwan form the parameterized and validated BRNM       #
#       * Draw conclusions and write the paper.                        #
#                                                                      #
# - I use similar parameters as used in the notebook                   #
# simple_metapopulation_SEIR_plus_mosquito_v0.0.1 in order to test     #
# real value parameters for chikungunya                                #
# - In this version I test the three patches model                     #
# - In this version I use the real conection data from openflights.org #
# to build up the mobility matrix 'M', I use real population sizes,    #
# of the Caribbena region of interest for the CHIKV outbreak           #
# - This is the same as the version 0.4 but a bit more clean           #
# - in this version I include a new node: Saint Martin and use the     #
# information given by the airport http://www.saintmartin-airport.com  #
# to update the conection mobility matrix. This should be important    #
# because the first reported local transmission in the Caribbean was   #
# reported precisely in Sanint Martin                                  #
# version 0.4.5:                                                       #
#   - add Tycho data for comparison and calibration of the theoretical #
#   - MODEL TWEAKER: this is a no-frills code that produces only the   #
#   theoretical aggragated curve and it compares against the (Tycho)   #
#   observed data. The goal is to have a clear nd simple code to       #
#   calibrate the theoretcial model with the observations by tweaking  #
#   the parameters                                                     #
#   - SCENARIOS:                                                       #
#       a) randomize the conections Saint Martin                       #
#          (first island infected) and loop for several realizations.  #
#          The goal is to check whether the pattern of sequnces is hold
#                                                                      #
#       b) Random seeding of the diseases in any island  while keeping
#          the topology and                                            #
#       several realizations - check whether the sequential            #
#       pattern holds                                                  #
#                                                                      #
#       c) connect Saint Martin with all of the other islands (Xavi's  #
#       idea -?)                                                       #
#                                                                      #
#       d) Random network with the same number of nodes and links of   #
#        the obeserved one.                                            #
#                                                                      #
#   version 0.6.1:                                                     #
#       a) Cleaning up existing code and refactoring funcions          #
#                                                                      #
########################################################################

########################################################################

__author__ = """carlos j dommar ~ carlos.dommar_at_isglobal.org"""
# Preambule:
import scipy.integrate as spi
import numpy as np
import pylab as pl
import pandas as pd
import networkx as nx
import seaborn
import time
__author__ = """carlos j dommar ~ carlos.dommar _at_isglobal.org"""


# Keep track of the tota ime of execution
start = time.time()

#############################################
# LOADING DATASETS AND BUILDING DATAFRAMES. #
#############################################
# all airports in the openflights.org dataset-
## use for the Lenovo laptop:
# airports = pd.read_csv('/home/carlos/projects/fpol/data/airports.csv')

## use for Ingold:
airports = pd.read_csv('./data/airports.csv')
# Load the outcome from the albert's script 'connection_matrix.py' (report path)
## use the following path for the Lenovo laptop:
#my_airports_df = \
#    pd.read_csv('/home/carlos/projects/chik-caribbean/non-in-repo/data00/island_connections_caribbean_3.csv')
## use the following pathe for Ingold:
my_airports_df = \
    pd.read_csv('~/projects/chik-caribbean/non-in-repo/data00/island_connections_caribbean_3.csv')
# cleaning off some unnecesary header off the dataframe
my_airports_df = my_airports_df.drop('Unnamed: 0',1)
my_airports_df.index = my_airports_df.columns

# Some cleaning and editing of the database is needed:
airports.ix[2821, 'Country'] = 'Bonaire'
airports.ix[2822, 'Country'] = 'Curacao'
airports.ix[2823, 'Country'] = 'Sint Eustatius'
airports.ix[2824, 'Country'] = 'Sint Maarten'
airports.ix[4137, 'Country'] = 'Saba'
airports.ix[5331, 'Country'] = 'Saint Barthelemy'

# Create a simple DF with airports with latitudes and longitudes:
columns = ['name', 'city', 'country', 'lat', 'lon']
index = my_airports_df.index
# DF for lat and long of selected airports:
latlonglobal_df = pd.DataFrame(index=index, columns=columns)
for i in latlonglobal_df.index:
    name = airports[airports['Airport ID'] == np.int(i)]['Name'].iloc[0]
    city = airports[airports['Airport ID'] == np.int(i)]['City'].iloc[0]
    country = airports[airports['Airport ID'] == np.int(i)]['Country'].iloc[0]
    lat = airports[airports['Airport ID'] == np.int(i)]['Latitude'].iloc[0]
    lon = airports[airports['Airport ID'] == np.int(i)]['Longitude'].iloc[0]
    latlonglobal_df.set_value(i, 'name', name)
    latlonglobal_df.set_value(i, 'city', city)
    latlonglobal_df.set_value(i, 'country', country)
    latlonglobal_df.set_value(i, 'lat', lat)
    latlonglobal_df.set_value(i, 'lon', lon)
del columns, index, lat, lon, i, name, city

print ('----------------------------------------------')
print ('The countries in the region of interest are the following:')
print ('')
print (latlonglobal_df['country'].drop_duplicates())
print ('')
print ('')
print ('They are ', latlonglobal_df['country'].drop_duplicates().shape[0], \
    ' countries')
# Add an attriute for 'political region' to the data frame:
latlonglobal_df['political_region'] = pd.Series(data=None)
# 'political_region' can e either. 'french', 'anglo', 'latin', 'ducth'
# print latlonglobal_df.head()


# now we assign the political regions
for i in latlonglobal_df['country']:
    if (i == 'Puerto Rico') or (i == 'Dominican Republic') \
            or (i =='Cuba') or (i == 'Colombia') or (i == 'Venezuela') \
            or (i == 'Mexico') or (i == 'Brazil') \
            or(i=='Nicaragua') or (i=='Panama') \
            or (i=='Costa Rica') or (i == 'Haiti'):
        latlonglobal_df.ix[latlonglobal_df['country'] \
                           == i, 'political_region'] = 'latin'

    elif (i == 'Jamaica') or (i == 'Cayman Islands') \
        or (i == 'Antigua and Barbuda') or (i == 'Grenada') \
        or (i == 'Barbados') or (i == 'Virgin Islands')  \
        or (i == 'Saint Kitts and Nevis') or (i == 'Anguilla') \
        or (i =='Trinidad and Tobago') or (i == 'Jamaica') \
        or (i=='British Virgin Islands')\
        or (i == 'Saint Vincent and the Grenadines') or (i =='Montserrat') \
        or (i=='Bahamas') or (i == 'Turks and Caicos Islands') \
        or (i =='Dominica') or (i =='Saint Lucia') or (i =='Montserrat') \
        or (i == 'Guyana'):
        latlonglobal_df.ix[latlonglobal_df['country'] \
                               == i, 'political_region'] = 'anglo'

    elif (i == 'Martinique') or (i =='Guadeloupe') or (i =='Saint Barthelemy') \
        or (i == 'French Guiana') or (i == 'Haiti'):
        latlonglobal_df.ix[latlonglobal_df['country'] \
                           == i, 'political_region'] = 'french'

    elif (i =='Aruba') or (i == 'Suriname') or (i == 'Bonaire') \
        or (i == 'Sint Eustatius') or (i == 'Curacao') \
        or (i == 'Sint Maarten') or (i == 'Saba'):
        latlonglobal_df.ix[latlonglobal_df['country'] \
                           == i, 'political_region'] = 'dutch'



##########################################################
##########################################################
# CONSTRUCTING A GRAPH FOR THE FLIGHT CONNECTION NETWORK #
##########################################################
##########################################################
# Remember to import networkx
# New method to create the graphs
G = nx.DiGraph(np.asarray(my_airports_df)) # direct method
new_nodes = np.asarray(my_airports_df.columns) # new nodes
mapping = {} # dict
count = 0
for i in G.node:
    mapping[i] = new_nodes[count]
    count = count + 1
del count, new_nodes
# new graph  with relabeld nodes:
G = nx.relabel_nodes(G, mapping)
# add some information to the nodes
for i in my_airports_df.index: # = my_aiports.columns
    G.add_node(i,
               airport_name = \
               latlonglobal_df[latlonglobal_df.index==i]['name'][0],

               city = \
               latlonglobal_df[latlonglobal_df.index==i]['city'][0],

               country = \
               latlonglobal_df[latlonglobal_df.index==i]['country'][0],

               lat = \
               latlonglobal_df[latlonglobal_df.index==i]['lat'][0],

               lon = \
               latlonglobal_df[latlonglobal_df.index==i]['lon'][0],

               political_region = \
               latlonglobal_df[latlonglobal_df.index==i]\
               ['political_region'][0],
               color = 'Gray',
               )
# remove isolates nodes:
G.remove_nodes_from(nx.isolates(G))

# setting the nodes some centrality attributes:
bb = nx.betweenness_centrality(G)
nx.set_node_attributes(G, 'betweeness', bb)
del bb

degree = nx.degree(G) # dict with nodes and corresponding degree

# print " "
# pl.figure(figsize=(5,5))
# nx.draw(G, node_color='gray', alpha=.60, node_size=200) #,
# pl.title('Graph of flight routes in The Caribbean region 3 (openflight.org)', \
#          fontsize=20)
# # pl.savefig("graph_region2.png",dpi=400)
# pl.show()

# Set color of nodes by political regions  for visualizations purposes
# Set the node attribute a different color per region:
for i in G.nodes():
    if G.node[i]['political_region'] == 'latin':
        G.node[i]['color'] = 'Peru'

    elif G.node[i]['political_region'] == 'dutch':
        G.node[i]['color'] = 'OrangeRed'

    elif G.node[i]['political_region'] == 'anglo':
        G.node[i]['color'] = 'DarkGreen'

    elif G.node[i]['political_region'] == 'french':
        G.node[i]['color'] = 'DarkBlue'

# Make a list with the color of the nodes per region using list comprehesion:
color_of_nodes = [G.node[i]['color'] for i in G.nodes()]



############################
# DATA GROUPED BY COUNTRY: #
############################
# make a list of the locations involved at the level of country/island:
countries = []
for i in G.nodes():
    countries.append(G.node[i].get('country'))

# drop duplicates countries:
# trick: use 'set'. A set is something that can't possibly have duplicates:
countries = list(set(countries))

# Correction: I need to add Saint Martin
countries.append('Saint Martin')


# The connecivity matrix is first computes as a Pandas
# Data Frame object ''
countries_df = pd.DataFrame(data=0, \
                            index=countries, columns=countries, dtype=int)

# Then
origin = str()
destiny = str()
for i in range(0, len(my_airports_df)):
    airport_ID = my_airports_df.index[i]
    if my_airports_df.iloc[i].sum() != 0:
        origin = latlonglobal_df[latlonglobal_df.index \
                                 == airport_ID]['country'][0]
        for j in range(0, len(my_airports_df)):
            if my_airports_df.iloc[i,j] != 0: #there's a connection in (i,j)
                destination_ID = my_airports_df.columns.values[j]
                destiny = latlonglobal_df[latlonglobal_df.index \
                                          == destination_ID]['country'][0]
                if origin != destiny:
                    #print 'Two different countries are connected!'
                    countries_df.ix[origin,destiny] \
                        = countries_df.ix[origin,destiny] \
                        + my_airports_df.iloc[i,j]
                destiny = list()
        origin = list()
# check that the diagonal of 'countries_df' must be zero:
for i in range(0,len(countries_df)):
    for j in range(0,len(countries_df)):
        if i == j and countries_df.iloc[i,j] != 0:
            print ('ALERT DIAGONAL NOT ZERO')
del i, j, origin, destiny, airport_ID


# Correcting Saint Martin connections:
# according the airport website http://www.saintmartin-airport.com
# Saint Martin has conection with Martinique, Saint Barthelemy,
# Guadalope, Saint Lucia, Santo Domingo (Puerto Rico), and Miami (USA):
# with Saint Barthelemy:  1
# with Martinique:        2
# with Saint Lucia:       1
# with Guadeloupe:        2
# with Puerto Rico:       1
# with Miami (USA):       2
countries_df.ix['Saint Martin', 'Saint Barthelemy'] = 1
countries_df.ix['Saint Barthelemy', 'Saint Martin'] = 1
countries_df.ix['Saint Martin', 'Martinique'] = 2
countries_df.ix['Martinique', 'Saint Martin'] = 2
countries_df.ix['Saint Martin', 'Saint Lucia'] = 1
countries_df.ix['Saint Lucia', 'Saint Martin'] = 1
countries_df.ix['Saint Martin', 'Guadeloupe'] = 2
countries_df.ix['Guadeloupe', 'Saint Martin'] = 2
countries_df.ix['Saint Martin', 'Puerto Rico'] = 1
countries_df.ix['Puerto Rico', 'Saint Martin'] = 1



# The connectivity matrix as a numpy array:
M = np.asarray(countries_df)
# New method to create the graphs
H = nx.DiGraph(M)  # direct method
mapping = {} # dict
count = 0
for i in range(0, len(countries)):
    mapping[i] = countries[count]
    count = count + 1
del count

## (but anyway I have done some coding before building the G network so I will re-use some code)
## G is the network by airport of the selected region and they seem to be 47
## on the other hand, H is the network by countries of the selected region and they are 30
H = nx.relabel_nodes(H, mapping)
# add attributes to H
for i in H.nodes():
    for j in G.nodes():
        if i == G.node[j].get('country'):
            # print 'yes, I found ', i, j
            H.add_node(i, color = G.node[j].get('color'), \
                       political_region = G.node[j].get('political_region'))
            break
#for Saint Martin:
H.add_node('Saint Martin', color = 'DarkBlue', political_region = 'french')

pop = pd.read_csv("./data/locations_populations.csv")
# Add population sizes and densities to the graph H:
for i in H.nodes():
    for j in pop['location']:
        if i == j:
            H.node[i]['pop_size'] = \
                np.float(pop[pop['location'] == j]['pop_size'])
            H.node[i]['density'] = \
                np.float(pop[pop['location'] == j]['density'])


# setting the nodes some centrality attributes:
bb = nx.betweenness_centrality(H)
nx.set_node_attributes(H, 'betweeness', bb)
del bb


# REMARK: There are airports from some countries that have no connection with \
# other airports in the chosen area. These may be airports in the borders of \
# the regions which connections goe sto some other locations outside the chosen\
# area and for that reason they appear as to be islotaed conountries. The \
# following code helps us to indentify these countries.
count = 0
for i in countries_df.index:
    # the country appears as having no connections with others in the area
    if countries_df.ix[i,:].sum() == 0:
        print ('The airports of the country  \
        of the chosen area seems to have no connection with the rest:')
        print(i)
        count = count + 1
print ('-----------------------------')
print ('There are ', count, ' such countries in the chosen area.')
print ('-----------------------------')
del count

# adding a tocuh of color:
color_of_nodes = [H.node[i]['color'] for i in H.nodes()]


# Simulation specs:
epi_classes = 5  # SEIRN
# patches = 3  # Number of patches or nodes in the SEIR-network
patches = nx.number_of_nodes(H)
leap = epi_classes * (patches - 1) + 1  # patches = is the number of patches
ND=365.0/1.      # number of days
ND=45.0*7        # number weeks times number of days
#TS=.10           # time step
TS=1.0           # time step



############################
# NETWORKED SEIR MODELLING #
############################
# Disease infection parameters:
# mu=1/(70*365.0)  # natality/mortality rate for a demographic model
beta=182.5/365.0 # transmission rate
beta = .78       # (modded)
sigma=1/4.032    # (IIP)
sigma=1/2.032    # (IIP) (modded)
gamma=1/7.194    # recovery rate
gamma=1/4.194    # recovery rate (modded)

beta = 1.4*beta
beta = np.round(beta, 2)
sigma = 1.2*sigma
#print ("sigma bfeore: ", sigma)
sigma = sigma  = np.round(sigma, 2)
#print ("sigma after: ", sigma)
gamma = 1.3*gamma
gamma = np.round(gamma, 2)

# Simulation specs:
epi_classes = 5  # SEIRN
# patches = 3  # Number of patches or nodes in the SEIR-network
patches = nx.number_of_nodes(H)
leap = epi_classes * (patches - 1) + 1  # patches = is the number of patches
ND=365.0/1.      # number of days
ND=45.0*7        # number weeks times number of days
#TS=.10           # time step
TS=1.0           # time step

# Mobility matrix:
# Convert aggregated by countries graph 'H' into a M numpy ndarray for speed:
M = nx.to_numpy_matrix(H)
dim = len(M)     # dimension of the connection matrix = number of nodes
N0 = np.zeros(shape=(dim,), dtype=float)
del dim

####################################!###########################################
################################################################################
# TEST: rescaling population sizes:                                            #
# This snippet calibrate the magnitude of the number of infectious cases       #
# by decreasing the population size exposed to the infection                   #
pop_rescale = .0012                                                            #
for i in H.nodes():                                                            #
    H.node[i]['pop_size'] = H.node[i]['pop_size'] * pop_rescale                #
beta = .69*beta                                                                #
beta = np.round(beta, decimals=2)
sigma = 1.3*sigma                                                              #
#print("sigma before 2: ", sigma)
sigma = np.round(sigma, decimals=2)
#print("sigma after 2:", sigma)
gamma = 1.4*gamma                                                              #
gamma = np.round(gamma, decimals=2)
m_factor = .56e-2*1                                                            #
M = M*m_factor # rescale matrix                                                #
################################################################################
####################################!###########################################

count = 0
for i in H.nodes():
    # N0[count] = H.node[i]['density'] # for densities
    N0[count] = H.node[i]['pop_size'] # for pop_sizes
    count = count + 1
del count

# initial vector
Y0 = np.zeros(epi_classes * patches)
# fill the initial vector:
count = 0
for n in N0:
    Y0[count] = n # for class S
    Y0[count+4] = n # for class N
    count = count + 5
del count

## correction for Saint Martin:
Y0[0] = Y0[0] - 1  # S = S - 1
Y0[2] = Y0[2] + 1  # I = I + 1
## #################################



#################################################################################
#################################################################################
# THE EQUATIONS:
def diff_eqs(Y0, t, M):
    # I need to keep track explicitely on the total population sizes
    S = Y0[0:leap:epi_classes]
    E = Y0[1:leap+1:epi_classes]
    I = Y0[2:leap+2:epi_classes]
    R = Y0[3:leap+3:epi_classes]
    N = Y0[4:leap+4:epi_classes]

    ## lines introduced by leo:
    #for i in range(0,len(I)):
    #            if I[i] > 1:


    dS = -beta*S*I/(S+E+I+R)
    dE = beta*S*I/(S+E+I+R) - sigma*E
    dI = sigma*E - gamma*I
    dR = gamma*I
    #print("sigma from inside the ode solver: ", sigma)

    # Update with migration rates:
    # For all dS's:
    for i in range(0,len(dS)):
        S_plus = (np.squeeze(np.asarray(M[:,i]))*S/(S+E+I+R)).sum()
        S_minus = M[i,:].sum()*S[i]/(S[i]+E[i]+I[i]+R[i])
        np.rint(S_plus)

        dS[i] = dS[i] \
                - S_minus \
                + S_plus
            #- M[i,:].sum()*S[i]/(S[i]+E[i]+I[i]+R[i])\
            #+ (np.squeeze(np.asarray(M[:,i]))*S/(S+E+I+R)).sum()

    # For all dE's:
    for i in range(0,len(dE)):
        E_plus = (np.squeeze(np.asarray(M[:,i]))*E/(S+E+I+R)).sum()
        E_minus = M[i,:].sum()*E[i]/(S[i]+E[i]+I[i]+R[i])\
        # threshold to simulate discretization:
        np.rint(E_plus)

        dE[i] = dE[i] \
                - E_minus \
                + E_plus
            #- M[i,:].sum()*E[i]/(S[i]+E[i]+I[i]+R[i])\
            #+ (np.squeeze(np.asarray(M[:,i]))*E/(S+E+I+R)).sum()

    # For all dI's:
    for i in range(0,len(dI)):
        I_plus = (np.squeeze(np.asarray(M[:,i]))*I/(S+E+I+R)).sum()
        I_minus = M[i,:].sum()*I[i]/(S[i]+E[i]+I[i]+R[i])
        np.rint(I_plus)

        dI[i] = dI[i] \
                - I_minus + I_plus
            #- M[i,:].sum()*I[i]/(S[i]+E[i]+I[i]+R[i])\
            #+ (np.squeeze(np.asarray(M[:,i]))*I/(S+E+I+R)).sum()

    # For all dR's
    for i in range(0,len(dR)):
        dR[i] = dR[i] \
            - M[i,:].sum()*R[i]/(S[i]+E[i]+I[i]+R[i])\
            + (np.squeeze(np.asarray(M[:,i]))*R/(S+E+I+R)).sum()

    f = [dS, dE, dI, dR, (dS+dE+dI+dR)]
    #f = np.rint([dS, dE, dI, dR, (dS+dE+dI+dR)])

    # put f in the right order:
    f_temp = list()
    for p in range(0, patches):
        for i in f:
            f_temp.append(i[p])
    f = f_temp # update f with the rightly sorted one
    f_temp = list() # reset f_temp
    return f # return to odeint

#################################################################################
#################################################################################



#################################################################################
#################################################################################
# THE EQUATIONS V.2.:
def diff_eqs_2(Y0, t, M):
    # I need to keep track explicitely on the total population sizes
    S = Y0[0:leap:epi_classes]
    E = Y0[1:leap+1:epi_classes]
    I = Y0[2:leap+2:epi_classes]
    R = Y0[3:leap+3:epi_classes]
    N = Y0[4:leap+4:epi_classes]
    ## lines introduced by leo:
    #for i in range(0,len(I)):
    #            if I[i] > 1:
    #print("sigma from inside the ode solver B: ", sigma)


    dS = -beta*S*I/(S+E+I+R)
    dE = beta*S*I/(S+E+I+R) - sigma*E
    dI = sigma*E - gamma*I
    dR = gamma*I

    # Update with migration rates:
    # For all dS's:
    for i in range(0,len(dS)):
        S_plus = (np.squeeze(np.asarray(M[:,i]))*S/(S+E+I+R)).sum()
        S_minus = M[i,:].sum()*S[i]/(S[i]+E[i]+I[i]+R[i])
        np.rint(S_plus)

        dS[i] = dS[i] \
                - S_minus\
                + S_plus
            #- M2[i,:].sum()*S[i]/(S[i]+E[i]+I[i]+R[i])\
            #+ (np.squeeze(np.asarray(M2[:,i]))*S/(S+E+I+R)).sum()

    # For all dE's:
    for i in range(0,len(dE)):
        E_plus = (np.squeeze(np.asarray(M[:,i]))*E/(S+E+I+R)).sum()
        E_minus = M[i,:].sum()*E[i]/(S[i]+E[i]+I[i]+R[i])\
        # threshold to simulate discretization:
        np.rint(E_plus)

        dE[i] = dE[i] \
                - E_minus \
                + E_plus
            #- M2[i,:].sum()*E[i]/(S[i]+E[i]+I[i]+R[i])\
            #+ (np.squeeze(np.asarray(M2[:,i]))*E/(S+E+I+R)).sum()

    # For all dI's:
    for i in range(0,len(dI)):
        I_plus = (np.squeeze(np.asarray(M[:,i]))*I/(S+E+I+R)).sum()
        I_minus = M[i,:].sum()*I[i]/(S[i]+E[i]+I[i]+R[i])
        np.rint(I_plus)

        dI[i] = dI[i] \
                - I_minus \
                + I_plus
            #- M2[i,:].sum()*I[i]/(S[i]+E[i]+I[i]+R[i])\
            #+ (np.squeeze(np.asarray(M2[:,i]))*I/(S+E+I+R)).sum()

    # For all dR's
    for i in range(0,len(dR)):
        R_plus = (np.squeeze(np.asarray(M[:,i]))*R/(S+E+I+R)).sum()
        R_minus = M[i,:].sum()*R[i]/(S[i]+E[i]+I[i]+R[i])
        np.rint(I_plus)
        dR[i] = dR[i] \
                - R_minus + R_plus
            #- M2[i,:].sum()*R[i]/(S[i]+E[i]+I[i]+R[i])\
            #+ (np.squeeze(np.asarray(M2[:,i]))*R/(S+E+I+R)).sum()

    f = [dS, dE, dI, dR, (dS+dE+dI+dR)]
    #f = np.rint([dS, dE, dI, dR, (dS+dE+dI+dR)])
    # Leonardo's suggestion:
    # It should wor because rint threshold to 1. values >= to .5 and to 0. values <.5
    #np.rint(f = [dS, dE, dI, dR, (dS+dE+dI+dR)])
    #f = np.rint([dS, dE, dI, dR, (dS+dE+dI+dR)])



    # put f in the right order:
    f_temp = list()
    for p in range(0, patches):
        for i in f:
            f_temp.append(i[p])
    f = f_temp # update f with the rightly sorted one
    f_temp = list() # reset f_temp
    return f # return to odeint

 
#################################################################################
#################################################################################


start = time.time()

# I run over a vector of differnte intervention dates (times)
intervention_dates = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]
intervention_dates = [1, 5, 10, 20, 40, 60, 80, 100, 120, 140, 150]
#intervention_dates = [1, 5, 40, 100]
#intervention_dates = [1, 40, 100]

#M = M*0.5 ## it seems it doesn not make  any effect here!!!
M2 = M * 0.0
interventions = [M, M*.9, M*.8, M*.7, M*.6,  M*.5, M*.4, M*.3, M*.2, M*.1, M*.0]

solution_list = [] # each element represents an intervention time from the list
                   # 'interventions'. Within each elemet there is a pd.Dataframe
                   # representing 10 intervention levels for that intervention
                   # date

for date in intervention_dates: # from 5, 10, 20, ... 150
    # SOLVE THE ODE:
    ## RESET INITIAL CONDITIONS TO ORIGINAL SET (from here)##
    ## Firstly I restart the values of the connectivity matrix to its originals:
    M = nx.to_numpy_matrix(H)
    m_factor = .56e-2*1
    M = M*m_factor # rescale matrix
    #"""

    # the intial condition vector Y0:
    count = 0
    for i in H.nodes():
        # N0[count] = H.node[i]['density'] # for densities
        N0[count] = H.node[i]['pop_size'] # for pop_sizes
        count = count + 1
    del count

    # initial vector
    Y0 = np.zeros(epi_classes * patches)
    # fill the initial vector:
    count = 0
    for n in N0:
        Y0[count] = n # for class S
        Y0[count+4] = n # for class N
        count = count + 5
    del count
    ## correction for Saint Martin:
    Y0[0] = Y0[0] - 1  # S = S - 1
    Y0[2] = Y0[2] + 1  # I = I + 1
    ## RESET INITIAL CONDITIONS TO ORIGINAL SET (up to here)##


    ND=45.0*7        # number weeks times number of days
    ############################
    ND_previous = date
    ############################

    #TS=.10           # time step
    TS=1.0           # time step
    t_start_previous = 0.0; t_end_previous = ND_previous; t_inc = TS
    t_range_previous = np.arange(t_start_previous, t_end_previous, t_inc)
    t_range = np.arange(0.0, ND, t_inc)

    #M=M*0.0
    # Here I run the previous epidemic before intervention
    #YSOL = spi.odeint(diff_eqs, Y0, t_range, args=(a2b,b2a))
    YSOL_0 = spi.odeint(diff_eqs, Y0, t_range_previous, (M,)) # YSOL_0 is the solution withou any modification
    YSOL_0 = YSOL_0*.15 # only 15% are observed or recorded


    ## Here I make a loop to run over different levels of decrease of flight traffic
    ## we cosnider this intervention scenarios as a GENERAL DECREASE of flight traffic
    ## by all locations t_start_post time after first case is reported


    t_star_post = t_end_previous; t_end_post = ND
    t_range_post = np.arange(t_end_previous, t_end_post, t_inc)

    # The new initial conditions are the conditions left by the connection matrix
    # is the new intervented one:
    Y0_post = YSOL_0[len(YSOL_0)-1, :]
    #Y0_post = Y0


    #del (diff_eqs)

    MAP = 'winter'
    #MAP = 'RdYlBu'
    cmap = pl.get_cmap(MAP)
    color = [cmap(i) for i in np.linspace(0, 1, len(interventions)+1)]
    pl.figure(figsize=(15,5))
    counter = 1
    # matrix solution, for each intervention date there are len(interventions)
    # levels:
    sol_matrix = np.zeros((len(t_range), len(interventions)))

    interv_num = 0
    for interv in interventions: # descending from 1.0,0.9,...,0.0
        YSOL_1 = spi.odeint(diff_eqs_2, Y0_post, t_range_post, (interv,)) # YSOL_0 is the solution without any modification
        #YSOL_1 = YSOL_1*.15 # only 15% are observed or recorded

        Total_Inf_0 = np.zeros(len(t_range_previous))
        count = 2
        for i in range(0,30):
            Total_Inf_0 = Total_Inf_0 + YSOL_0[:,count]
            count = count + 5
        del  count, i

        Total_Inf_1 = np.zeros(len(t_range_post))
        count = 2
        for i in range(0,30):
            Total_Inf_1 = Total_Inf_1 + YSOL_1[:,count]
            count = count + 5
        del  count,i

        # I need to drop the first row of Y_SOL1 since is redundat with the last row of YSOL_0
        YSOL_Total = np.concatenate((YSOL_0, YSOL_1[1:len(YSOL_1),:]))
        YSOL_Total = np.concatenate((YSOL_0, YSOL_1))
        Total_Inf_total = np.zeros(len(np.concatenate((t_range_previous, t_range_post))))
        count = 2
        for i in range(0,30):
            Total_Inf_total = Total_Inf_total + YSOL_Total[:,count]
            count = count + 5
        del  count, i

        # store the solution by intervention level:
        sol_matrix[:, interv_num] = Total_Inf_total
        print("interv_num =", interv_num)
        print("shape of 'Total_Inf_total' =", Total_Inf_total.shape)


        # Plotting both total PRE and POST INTERVENTION with different colors
        pl.title('Model curve of PRE (red) + POST (color gradient) INTERVENTION -aggregated infecteds- day %d' %ND_previous, fontsize=20)
        pl.plot(t_range_previous, Total_Inf_0, linewidth=7, alpha=1.0, color="Crimson")
        #if interv.all() == M.all():
        #    pl.plot(t_range_post, Total_Inf_1, linewidth=4, alpha=.7, color="black")
        pl.plot(t_range_post, Total_Inf_1, linewidth=4, alpha=.7, color=color[counter])
        pl.xlabel('time (days)', fontsize=20)
        pl.legend(loc=1)

        counter = counter + 1

        interv_num = interv_num + 1


    solution_list.append(sol_matrix)
    del counter

    ### PLotting
    ##pl.legend(['M', 'M*.9', 'M*.8', 'M*.7', 'M*.6',  'M*.5', 'M*.4', 'M*.3', 'M*.2', 'M*.1', 'M*.0'], loc='upper left')
    #pl.plot(ND_previous, Total_Inf_0[ND_previous*10 - 1] + 40, 'v', color='black', markersize=15) # marking the day of intervention
    ##pl.savefig("aggregated_infecteds_PRE_plus_POST_INTERVENTION_day_(%d).png" %ND_previous, dpi=400)
    ##pl.savefig("aggregated_infecteds_PRE_plus_POST_INTERVENTION_day_(%d).jpg" %ND_previous, pi=400)
    #pl.show(block=False)

################################################
## save and retrieve the solution list:
# save the solution list:
import pickle
with open("my_solution_list.txt", "wb") as fp:
    pickle.dump(solution_list, fp)

## to read the saved list:
#with open("my_solution_list.txt", "rb") as fp:
#    solution_list = pickle.load(fp)
################################################


## Plotting summary plot with intervention results

count_list = 0

# vector to store peaks
peak_vec = np.zeros(solution_list[0].shape[1])

# vector to store days to peaks
days2peak_vec = np.zeros(solution_list[0].shape[1])

# vector to store accumulated infecteds
acc_infec_vec = np.zeros(solution_list[0].shape[1])

pl.figure(figsize=(12,6))
for date in intervention_dates:
    print("DATE = ", date)

    peak_vec = np.asarray( [max(solution_list[count_list][:,i])
                            for i in np.arange(0, solution_list[count_list].shape[1])] )
    print("peak_vec = ", peak_vec)
    pl.subplot(3,3,1)
    pl.scatter(np.repeat(date, solution_list[0].shape[1]), peak_vec, color=color)
    pl.ylabel("Infecteds at largest peak")
    pl.xlabel("Intervention day")


    days2peak_vec = np.asarray( [ np.argmax(solution_list[count_list][:,i])
                            for i in np.arange(0, solution_list[count_list].shape[1])] )
    days2peak_vec = days2peak_vec/10.
    print("days2preak = ", days2peak_vec )
    pl.subplot(3,3,5)
    pl.scatter(np.repeat(date, solution_list[0].shape[1]), days2peak_vec, color=color)
    pl.ylabel("Days to largest peak")
    pl.xlabel("Intervention day")



    acc_infec_vec = np.asarray([ sum( solution_list[count_list][0:np.argmax(solution_list[count_list][:,i]),i] )
                                 for i in np.arange(0, solution_list[count_list].shape[1] )])
    print("acc_infect_vect = ", acc_infec_vec)
    print("")
    pl.subplot(3,3,9)
    pl.scatter(np.repeat(date, solution_list[0].shape[1]), acc_infec_vec, color=color)
    pl.ylabel("Accum. infect.")
    pl.xlabel("Intervention day")

    pl.tight_layout()
    pl.show(block=False)
    count_list = count_list + 1

#pl.savefig("figures/intervention_summary.png",dpi=400)
#pl.savefig("figures/intervention_summary.pdf")



end = time.time()
print ("----------")
print ("Intervention made at day: ", ND_previous)
print ("execution time: ", end - start, " secs")
print ("converted to minutes: ", (end - start)/60, "minutes")
print ("----------")
del start, end
