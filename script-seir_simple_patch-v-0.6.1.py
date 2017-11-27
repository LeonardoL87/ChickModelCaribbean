#!/usr/bin/env python
# SEIR without demography (no vital dynamics)
# VERSION NOTEBOOK
##################################################################
#    This  porgrame started as the PYTHON version of program 2.6 #
# from page 41 of                                                #
# "Modeling Infectious Disease in humans and animals"            #
# by Keeling & Rohani                                            #
#	                                							 #
# It is the SEIR epidemic with no births nor deaths.             #l
##################################################################
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
#   version 0.6.1:
#       a) Cleaning up existing code and refactoring funcions
#                                                                      #
########################################################################

########################################################################

__author__ = """carlos j dommar ~ carlos.dommar_at_ ic3.cat"""
# Preambule:
import os
import scipy.integrate as spi
import numpy as np
import pylab as pl
import pandas as pd
import networkx as nx
import seaborn
import time
#%matplotlib inline

start = time.time()
###############################################################################
########################## MY FUNCTIONS #######################################
###############################################################################
###############################################################################
###############################################################################


###############################################################################
# F_AIRPORTS(DATA_PROJECT_FOLDER)
def f_airports(data_project_folder):
   """
   this function receives a path variable where the data is located
   and gives back a pandas dataframe with airports by label
   and their connections.
   """
   ## use for both Ingold and Lenovo:
   airports = pd.read_csv('./data/airports.csv')

   # Load the outcome from the albert's script 'connection_matrix.py' (report path)
   ## use the following path for the Lenovo laptop:
   airports_df = pd.read_csv('./data/island_connections_caribbean_3.csv')
   # cleaning off some unnecesary header off the dataframe
   airports_df = airports_df.drop('Unnamed: 0',1)
   airports_df.index = airports_df.columns

   return airports_df


#################################################################################
# F_LATLON_DF()
def f_latlon_df():
    """
    this function recieves nothing (the path to the data is identical both in
    Ingold and in my laptop), and gives back a pandas dataframe  with the
    countries, their number tag from the original data frame, ther coordinates,
    political, region, etc.
    """
    ## use for both Ingold and Lenovo:
    airports = pd.read_csv('./data/airports.csv')
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
    print ('They are '), latlonglobal_df['country'].drop_duplicates().shape[0], \
        (' countries')
    # Add an attriute for 'political region' to the data frame:
    latlonglobal_df['political_region'] = pd.Series(data=None)
    # 'political_region' can e either. 'french', 'anglo', 'latin', 'dutch'
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

    return latlonglobal_df


#################################################################################
# F_GRAPH(DF, LATLONGLOBAL_DF)
def f_graph(df, latlonglobal_df): # rough network model
    """
#    this  function receives 'df', a pandas' dataframe  coming out the function
#    'my_aiports()', and a pandas's dataframe coming out the function
#    'latlong_df()' and gives back  a networkx object with teh final network model
#    for the flight connection of the area selected in the Albert's script
#    the funtion returns a networkx  object with the model of the observed
#    flight connection

    """
    # New method to create the graphs
    G = nx.DiGraph(np.asarray(df)) # direct method
    new_nodes = np.asarray(df.columns) # new nodes
    mapping = {} # dict
    count = 0
    for i in G.node:
        mapping[i] = new_nodes[count]
        count = count + 1
    del count, new_nodes
    # new graph  with relabeld nodes:
    G = nx.relabel_nodes(G, mapping)
    # add some information to the nodes
    for i in df.index: # = my_aiports.columns
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
    return G



#################################################################################
# F_DATAFRAME(MY_G, MY_AIRPORTS_DF, MY_LATLONG_DF)
def f_dataframe(my_g, my_airports_df, my_latlong_df):
    """
    this function receives:
        - a networkx object (usually coming out the function 'f_graph()')
        - a pandas data frame object with all selected airports (usually comming
        out the function 'f_airports()')
        - a pandas data frame object with all selected airports (usually coming
        out of the function 'f_latlon_df()')
    this function gives back: a pandas data frame object with corrected and
    aggregated connections by countries rather than by airports
    """
    ############################
    # DATA GROUPED BY COUNTRY: #
    ############################
    # make a list of the locations involved at the level of country/island:
    countries = []
    for i in my_g.nodes():
        countries.append(my_g.node[i].get('country'))

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
            origin = my_latlong_df[my_latlong_df.index \
                                     == airport_ID]['country'][0]
            for j in range(0, len(my_airports_df)):
                if my_airports_df.iloc[i,j] != 0: #there's a connection in (i,j)
                    destination_ID = my_airports_df.columns.values[j]
                    destiny = my_latlong_df[my_latlong_df.index \
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
    #return  countries_df


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

    return countries_df


################################################################################
# F_REDUCED_G(MYDF, MY_G)
def f_reduced_g(mydf, my_g, pop_rescale):
    """
    this function receives in the following order:
        - mydf: pandas data frame with the definite countries and their
        connections
        - my_g: netwokx object with connection of ALL airports
    this function gives back:
        - a networkx object with the definite (reduced) network model by country.
    """
    #mydf = countries_df  # finally the mydf dataframe with the selected area
    # The connectivity matrix as a numpy array:
    M = np.asarray(mydf)
    # new more clear line --- I had to transpose the matrix to get H right:
    H = nx.from_numpy_matrix(np.matrix.transpose(M))
    mapping = {} # dict
    count = 0

    # I need this piece of code again:
    countries = []
    for i in my_g.nodes():
        countries.append(my_g.node[i].get('country'))
    # drop duplicates countries:
    # trick: use 'set'. A set is something that can't possibly have duplicates:
    countries = list(set(countries))
    # Correction: I need to add Saint Martin
    countries.append('Saint Martin')

    for i in range(0, len(countries)):
        mapping[i] = countries[count]
        count = count + 1
    del count

    H = nx.relabel_nodes(H, mapping)
    # add attributes to H
    for i in H.nodes():
        for j in my_g.nodes():
            if i == my_g.node[j].get('country'):
                # print 'yes, I found ', i, j
                H.add_node(i, color = my_g.node[j].get('color'), \
                           political_region = my_g.node[j].get('political_region'))
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

    # TEST: rescaling population sizes:                                         #
    # This snippet calibrate the magnitude of the number of infectious cases    #
    # by decreasing the population size exposed to the infection                #
    for i in H.nodes():                                                         #
        H.node[i]['pop_size'] = H.node[i]['pop_size'] * pop_rescale


    # setting the nodes some centrality attributes:
    bb = nx.betweenness_centrality(H)
    nx.set_node_attributes(H, 'betweeness', bb)
    del bb, i, j

    # REMARK: There are airports from some countries that have no connection with \
    # other airports in the chosen area. These may be airports in the borders of \
    # the regions which connections goe sto some other locations outside the chosen\
    # area and for that reason they appear as to be islotaed conountries. The \
    # following code helps us to indentify these countries.
    count = 0
    for i in mydf.index:
        # the country appears as having no connections with others in the area
        if mydf.ix[i,:].sum() == 0:
            print ('The airports of the country '), i, \
                (' of the chosen area seems to have no connection with the rest.')
            count = count + 1
    print ('-----------------------------')
    print ('There are '), count, (' such countries in the chosen area.')
    print ('-----------------------------')
    del count, i

    # adding a tocuh of color:
    color_of_nodes = [H.node[i]['color'] for i in H.nodes()]

    return H


#################################################################################
# THE SEIR EQUATIONS: F_DIFF_EQS(YO, T, M, BETA, GAMMA, SIGMA)
def f_diff_eqs(Y0, t, epi_classes, M, beta, sigma, gamma):
    patches = M.shape[0]
    leap = epi_classes * (patches - 1) + 1  # patches = is the number of patches
    # I need to keep track explicitely on the total population sizes
    S = Y0[0:leap:epi_classes]
    E = Y0[1:leap+1:epi_classes]
    I = Y0[2:leap+2:epi_classes]
    R = Y0[3:leap+3:epi_classes]
    N = Y0[4:leap+4:epi_classes]  # never used!

    dS = -beta*S*I/(S+E+I+R)
    dE = beta*S*I/(S+E+I+R) - sigma*E
    dI = sigma*E - gamma*I
    dR = gamma*I

    # Update with migration rates:

    # For all dS's:
    for i in range(0,len(dS)):
        S_plus = (np.squeeze(np.asarray(M[:,i]))*S/(S+E+I+R)).sum()
        S_minus = M[i,:].sum()*S[i]/(S[i]+E[i]+I[i]+R[i])\

        dS[i] = dS[i] \
            #- M[i,:].sum()*S[i]/(S[i]+E[i]+I[i]+R[i])\
            #+ (np.squeeze(np.asarray(M[:,i]))*S/(S+E+I+R)).sum()
        - S_minus \
            + S_plus


    # For all dE's:
    for i in range(0,len(dE)):
        E_plus = (np.squeeze(np.asarray(M[:,i]))*E/(S+E+I+R)).sum()
        E_minus = M[i,:].sum()*E[i]/(S[i]+E[i]+I[i]+R[i])\

        ## threshold to simulate discretization:
        #if E_plus < 1E-5:
        #    E_plus = 0.0
        ###E_plus = 0.0
        #print ""
        #print ""
        #print ""
        #print ("dE = ",E_plus)
        #print ""
        #print ""

        dE[i] = dE[i] \
            #- M[i,:].sum()*E[i]/(S[i]+E[i]+I[i]+R[i])\
            #+ (np.squeeze(np.asarray(M[:,i]))*E/(S+E+I+R)).sum()
        - E_minus \
            + E_plus


    # For all dI's:
    for i in range(0,len(dI)):
        I_plus = (np.squeeze(np.asarray(M[:,i]))*I/(S+E+I+R)).sum()
        I_minus = M[i,:].sum()*I[i]/(S[i]+E[i]+I[i]+R[i])


        ## threshold to simulate discretization:
        #if I_plus < 1E-9:
        #    I_plus = 0.0
        #I_plus = 0.0
        ##print ""
        ##print ""
        ##print ""
        ##print ("dI = ",I_plus)
        ##print ""
        ##print ""

        dI[i] = dI[i] - I_minus + I_plus
            #- M[i,:].sum()*I[i]/(S[i]+E[i]+I[i]+R[i])\
            #+ (np.squeeze(np.asarray(M[:,i]))*I/(S+E+I+R)).sum()
            


    # For all dR's
    for i in range(0,len(dR)):
        dR[i] = dR[i] \
            - M[i,:].sum()*R[i]/(S[i]+E[i]+I[i]+R[i])\
            + (np.squeeze(np.asarray(M[:,i]))*R/(S+E+I+R)).sum()


    #f = [dS, dE, dI, dR, (dS+dE+dI+dR)]
    np.rint(f = [dS, dE, dI, dR, (dS+dE+dI+dR)])
    ### NUEVA VERSION PRUEBA ############
    ### NUEVA VERSION PRUEBA ############
    ### NUEVA VERSION PRUEBA ############
    ### NUEVA VERSION PRUEBA ############
    # put f in the right order:
    f_temp = list()
    patches = M.shape[0]  #number of patches
    for p in range(0, patches):
        for i in f:
            f_temp.append(i[p])
    f = f_temp # update f with the rightly sorted one
    f_temp = list() # reset f_temp

    return f # return to odeint


#################################################################################
# THE SEIR EQUATIONS (BY NODES): ODE_EQ(Y0_NODE, T, BETA, SIGMA, GAMMA)
def ode_eq(Y0_node, t, beta, sigma, gamma):
    # I need to keep track explicitely on the total population sizes
    S = Y0_node[0]
    E = Y0_node[1]
    I = Y0_node[2]
    R = Y0_node[3]

    dS = -beta*S*I/(S+E+I+R)
    dE = beta*S*I/(S+E+I+R) - sigma*E
    dI = sigma*E - gamma*I
    dR = gamma*I


    f = [dS, dE, dI, dR, (dS+dE+dI+dR)]

    return f # return to odeint
 

#################################################################################
# F_PLOT_OBSERVED_CHIK_INDV(DATA_PROJECT_FOLDER)
def f_plot_observed_chik_indv(data_project_folder):
    """
    this function receives:
        - a path variable 'data_project_folder',
    and it gives back:
        - a pandas.core.frame.DataFrame object that works as a dictionary where 
        the keys are the name of the nodes/locations and values as time series
        objects  with oserved chik cases
    """
    df = pd.read_csv('./data/Chikungunya_file_28August2014_Tycho_b.csv')

    # reaplce the cells with '\N' by 0:
    df.replace(r'\N', 0, inplace=True) # http://stackoverflow.com/questions/1347791/unicode-error-unicodeescape-codec-cant-decode-bytes-cannot-open-text-file

    # Getting the types right:
    df[['YEAR', 'WEEK', 'Week_conf_cases']] = df[['YEAR', 'WEEK', 'Week_conf_cases']].astype(int)
    df['Cum_sus_cases'] = df['Cum_sus_cases'].astype(int)
    df['Week_conf_cases'] = df['Week_conf_cases'].astype(int)
    df['Cum_conf_cases'] = df['Cum_conf_cases'].astype(int)
    df['Week_import_cases'] = df['Week_import_cases'].astype(int)
    df['Cum_import_cases'] = df['Cum_import_cases'].astype(int)
    df['Week_deaths'] = df['Week_deaths'].astype(int)
    df['Cum_deaths'] = df['Cum_deaths'].astype(int)
    df['Population'] =  df['Population'].astype(int)
    df['Week_Sus_Cases'] = df['Week_Sus_Cases'].astype(int)
    df['COUNTRY'] = df['COUNTRY'].astype(str)
    df['Incidence_sus_cases'] = df['Incidence_sus_cases'].astype(float)
    df['Incidence_conf_cases'] = df['Incidence_conf_cases'].astype(float)

    ## transforming negative data to absolute value:
    ## Example:
    #df.loc[df['COUNTRY'] == 'Virgin Islands, U.S.', 'Week_conf_cases']
    #virgin_islands_us = df.loc[df['COUNTRY'] == 'Virgin Islands, U.S.', 'Week_conf_cases']
    # function (inside a function! :-o):
    def neg2abs(x):
        if x < 0:
            return np.abs(x)
        else:
            return x

    #virgin_islands_us.apply(neg2abs)
    ## THERE'S TWO OPTIONS:
    ## 1. take the absolute value of all negative\
    ##    observations (assuming the the negative symbol is a mistake)
    df.ix[df.Week_conf_cases < 0, 'Week_conf_cases'] = \
        df.ix[df.Week_conf_cases < 0, 'Week_conf_cases'].apply(neg2abs)

    ## 2. set all negative observations to zero (assuming it is a mistake):
    #df.ix[df.Week_conf_cases < 0, 'Week_conf_cases'] = 0


    ##################################################
    ## TRANSFORMING ISO WEEK NUMBERS INTO PANDAS TIME:
    ##################################################
    from isoweek import Week
    #test:
    d = Week(2011, 51).sunday()
    #print d
    #print 'TRANFORMING ISO WEEK NUMBERS INTO PANDAS TIME: test OK'
    del d

    # CREATING PROPPER TIME SERIES OF OBSERVED DATA FROM THE TYCHO DATA BASE:
    ## create two lists: 'country' with country names and 'all_dates' with corresponding
    ## time indexes for observations
    count =  0
    country = []
    all_dates = []
    for x in df['COUNTRY'].drop_duplicates():
        if x != 0: # horrible hack here
            iso_week = np.asarray(df[df['COUNTRY'] == x][['YEAR','WEEK']])
            fecha = [Week(i[0], i[1]).sunday() for i in iso_week]
            # print count, ' -- ', x, fecha
            country.append(x)
            all_dates.append(fecha)
            count = count + 1
            # print ''
            del fecha
    del iso_week

    # create dictionary to relate country and time indexes
    dates_dict = {}
    countries_and_dates = zip(country, all_dates)
    for country, all_dates in countries_and_dates:
        dates_dict[country] = all_dates
    del countries_and_dates, all_dates


    # now the same but for a dictionary with the time series objects:
    count = 0
    country = []
    time_series = []
    for x in df['COUNTRY'].drop_duplicates():
        if x != 0: # horrible hack
            ts = pd.Series(np.asarray(df[df['COUNTRY'] == \
                                x]['Week_conf_cases']), dates_dict[x])
            # print count, ts
            country.append(x)
            time_series.append(ts)
            count = count + 1
            del ts

    # create dictionary to relate country and time series
    time_series_dict = {}
    countries_and_ts = zip(country, time_series)
    for country, time_series in countries_and_ts:
        time_series_dict[country] = time_series
    del countries_and_ts, time_series

    ## PROPERLY PLOT.SUM OF TIME SERIES WITH DIFFERENT LENGTHS
    # if the series are in the 'time_series_dict' dictionary:
    frame = pd.DataFrame(time_series_dict)
    # print frame.head()

    # Since different length are, in this particular, the result of
    # different starting times of the series, the appearing NaNs
    # mean that no cases where recorded in that dates. So we make them equals to zero:
    frame.fillna(0, inplace=True)
    # print frame.head()
    return(frame)


#################################################################################
# F_AGGREGATED_TS_OBSERVED_CHIK(DATA_PROJECT_FOLDER)
def f_aggregated_ts_observed_chik(data_project_folder):
    """
    this function receives a path variable and gives back a ts object (pandas
    time series) with the aggregated incidence data (first stages)
    reported by Tycho for the 2013 Caribbean outbreak of CHIKV
    """
    df = pd.read_csv('./data/Chikungunya_file_28August2014_Tycho_b.csv')

    # reaplce the cells with '\N' by 0:
    df.replace(r'\N', 0, inplace=True)

    # Getting the types right:
    df[['YEAR', 'WEEK', 'Week_conf_cases']] = df[['YEAR', 'WEEK', 'Week_conf_cases']].astype(int)
    df['Cum_sus_cases'] = df['Cum_sus_cases'].astype(int)
    df['Week_conf_cases'] = df['Week_conf_cases'].astype(int)
    df['Cum_conf_cases'] = df['Cum_conf_cases'].astype(int)
    df['Week_import_cases'] = df['Week_import_cases'].astype(int)
    df['Cum_import_cases'] = df['Cum_import_cases'].astype(int)
    df['Week_deaths'] = df['Week_deaths'].astype(int)
    df['Cum_deaths'] = df['Cum_deaths'].astype(int)
    df['Population'] =  df['Population'].astype(int)
    df['Week_Sus_Cases'] = df['Week_Sus_Cases'].astype(int)
    df['COUNTRY'] = df['COUNTRY'].astype(str)
    df['Incidence_sus_cases'] = df['Incidence_sus_cases'].astype(float)
    df['Incidence_conf_cases'] = df['Incidence_conf_cases'].astype(float)

    ## transforming negative data to absolute value:
    ## Example:
    #df.loc[df['COUNTRY'] == 'Virgin Islands, U.S.', 'Week_conf_cases']
    #virgin_islands_us = df.loc[df['COUNTRY'] == 'Virgin Islands, U.S.', 'Week_conf_cases']
    # function (inside a function! :-o):
    def neg2abs(x):
        if x < 0:
            return np.abs(x)
        else:
            return x

    #virgin_islands_us.apply(neg2abs)
    ## THERE'S TWO OPTIONS:
    ## 1. take the absolute value of all negative\
    ##    observations (assuming the the negative symbol is a mistake)
    df.ix[df.Week_conf_cases < 0, 'Week_conf_cases'] = \
        df.ix[df.Week_conf_cases < 0, 'Week_conf_cases'].apply(neg2abs)

    ## 2. set all negative observations to zero (assuming it is a mistake):
    #df.ix[df.Week_conf_cases < 0, 'Week_conf_cases'] = 0

    ##################################################
    ## TRASFORMING ISO WEEK NUMBERS INTO PANDAS TIME:
    ##################################################
    from isoweek import Week
    #test:
    d = Week(2011, 51).sunday()
    #print d
    #print 'TRANSFORMING ISO WEEK NUMBERS INTO PANDAS TIME: test OK'
    del d

    # CREATING PROPPER TIME SERIES OF OBSERVED DATA FROM THE TYCHO DATA BASE:
    ## create two lists: 'country' with country names and 'all_dates' with corresponding
    ## time indexes for observations
    count =  0
    country = []
    all_dates = []
    for x in df['COUNTRY'].drop_duplicates():
        if x != 0: # horrible hack here
            iso_week = np.asarray(df[df['COUNTRY'] == x][['YEAR','WEEK']])
            fecha = [Week(i[0], i[1]).sunday() for i in iso_week]
            # print count, ' -- ', x, fecha
            country.append(x)
            all_dates.append(fecha)
            count = count + 1
            # print ''
            del fecha
    del iso_week

    # create dictionary to relate country and time indexes
    dates_dict = {}
    countries_and_dates = zip(country, all_dates)
    for country, all_dates in countries_and_dates:
        dates_dict[country] = all_dates
    del countries_and_dates, all_dates

    # now the same but for a dictionary with the time series objects:
    count = 0
    country = []
    time_series = []
    for x in df['COUNTRY'].drop_duplicates():
        if x != 0: # horrible hack
            ts = pd.Series(np.asarray(df[df['COUNTRY'] == \
                                x]['Week_conf_cases']), dates_dict[x])
            # print count, ts
            country.append(x)
            time_series.append(ts)
            count = count + 1
            del ts

    # create distionary to relate country and time series
    time_series_dict = {}
    countries_and_ts = zip(country, time_series)
    for country, time_series in countries_and_ts:
        time_series_dict[country] = time_series
    del countries_and_ts, time_series

    ## PROPERLY PLOT.SUM OF TIME SERIES WITH DIFFERENT LENGTHS
    # if the series are in the 'time_series_dict' dictionary:
    frame = pd.DataFrame(time_series_dict)

    # Since different length are, in this particular, the result of
    # different starting times of the series, the appearing NaNs
    # mean that no cases where recorded in that dates. So we make them equals to zero:
    frame.fillna(0, inplace=True)
    # print frame.head()

    ## Print test:
    ## [frame[i].plot(linewidth=4, alpha=.35, \
    ##                figsize=(15,7)) for i in df['COUNTRY'].drop_duplicates()]
    #
    # Aggregate all oberved (Tycho) time series
    initial_country = df['COUNTRY'].drop_duplicates()[0]
    TS_total = frame[initial_country]

    for i in df['COUNTRY'].drop_duplicates():
        if i != initial_country:
            TS_total = TS_total + frame[i]

    del country
    return (TS_total)


#################################################################################
# F_YO_CONSTRUCTOR(H node)
def f_yo_constructor(H, node):
    """
    this function receives:
        - a netowrkx object that represents the oberved connection flights and
        the name (string) of the node where to start the infection
    this functiongives back:
        - an Y0 array with the initial condition for ode solver.
    """
    # initial vector
    patches = H.number_of_nodes()*epi_classes # # of patches
    Y0 = np.zeros(patches)
    # fill the initial vector:
    N0 = np.zeros(shape=(H.number_of_nodes()), dtype=float)

    count = 0
    for i in H.nodes():
        # N0[count] = H.node[i]['density'] # for densities
        N0[count] = H.node[i]['pop_size'] # for pop_sizes
        count = count + 1
    del count


    count = 0
    for n in N0:
        Y0[count] = n # for class S
        Y0[count+4] = n # for class N
        count = count + 5
    del count, N0

    ## correction for Saint Martin:
    # vector position for the asked node:
    pos = H.nodes().index(node)
    s_pos = epi_classes * pos
    i_pos = s_pos + 2

    Y0[s_pos] = Y0[s_pos] - 1  # S = S - 1
    Y0[i_pos] = Y0[i_pos] + 1  # I = I + 1
    return(Y0)


#################################################################################
## def f_sample_solution(my_observed_chik_ts, my_sol)
def f_sample_solution(my_observed_chik_ts, my_sol):
    """
    this function recieves:
        - a pandas time series object, representing the observed data
        - a pandas time series object, representing the solution
    this function gives back:
        - a pandas time series with a reduced vector solution with only data
        points corresponding to the obeserved datae for comparison purposes
    """
    ## need to transform the ts index to a DatetimeIndex to be able to extract
    ## later the values by date
    #my_observed_chik_ts.index = pd.DatetimeIndex(my_observed_chik_ts.index)
    my_reduced_sol = pd.Series()
    for i in my_observed_chik_ts.index:
        my_reduced_sol[i] = my_sol[str(i)].mean()

    return(my_reduced_sol)


#################################################################################
## def f_vector_distance(vector1, vector2)
def f_vector_distance(vector1, vector2):
    """
    this function compute the distance between any two real-valued vectors
    """
    dist = 0
    #while len(vector1) == len(vector2):
    for i in np.arange(0, len(vector1)):
        dist = dist + (vector1[i] - vector2[i])**2

    #else:
    #    print ("vector 1 and vector tw of the f_vector_distance functionmust be of equal length")

    return(np.sqrt(dist))


#################################################################################
# r_squared (y_obs, y_model)
def r_squared (y_obs, y_model):  # according to wikipedia
    """
    NOTE: seemingly this function should not be applied to nonlinear models (!)
    this function receives:
        - an observed vector, y_obs
        - a modeled or estimated vector, y_model.
    this funtions gives back the Coeffiecient of Determination  according the
    Wikipedia definition
    """
    SSTot = 0. # total sum of squares
    SSReg = 0. # regression sum of squares
    SSRes = 0. # regression sum of squares

    y_mean = y_obs.mean() # mean of the observed data

    for i in  np.arange(0, len(y_obs)):
        SSTot = SSTot +  (y_obs[i] - y_mean)**2
    del(i)
    #print "SSTot=", SSTot

    for i in np.arange(0, len(y_model)):
        SSReg = SSReg + (y_model[i] - y_mean)**2
    del(i)
    #print "SSReg=", SSReg

    for i in np.arange(0, len(y_obs)):
        SSRes = SSRes + (y_obs[i] - y_model[i])**2
    del(i)
    #print "SSRes=", SSRes

    R_square = 1 - SSRes/SSTot

    return(R_square)


#################################################################################
# MINIMUM_SQUARE_ERROR(MY_OBSERVED_CHIK_TS, H0, Y0, T_START, T_END, T_INC,
# T_RANGE,EP_CLASSES, BETA, SIGMA, GAMMA)
def minimum_square_error(my_observed_chik_ts,
                         H, Y0, t_start, t_end, t_inc, t_range, epi_classes,
                         init_beta, init_sigma, init_gamma,
                         init_bound, param_partition, recur, experiment):
    """
    This functions takes:
        - my_observed_chik_ts, a pandas time series object representing the
        aggregated beserved data of chik incidencein the Caribea outbreak
        reported by Tycho,
        - H, a networkx object that represents the oberved fliht connection for
        the area on those yeatas (2013) according openflight.org,
        - Y0, t_start, t_end, t_inc, t_range, conditions needed for teh ODE
        solver,
        - beta, sigma, and gamma , teh epi parameters to optimze by minimum mean
        squares
    This function gives back:
        - a tuple consisting of the  parameters (beta, sigma, gamma) that
        produces the solution with the shortet distance to the observed data
        for the corresponding discretization  of the theoretical solution
    """

    #start = time.time()

    a = init_beta - init_bound
    b = init_beta + init_bound
    if a < 0:
        a = 0
    beta_space = np.linspace(a, b, param_partition)
    #print "Beta_space = ", beta_space
    del (a,b)
    #-----------------

    a = init_sigma - init_bound
    b = init_sigma + init_bound
    if a < 0:
        a = 0
    sigma_space = np.linspace(a, b, param_partition)
    #print "Sigma_space = ", sigma_space
    del(a,b)
    #-----------------

    a = init_gamma - init_bound
    b = init_gamma + init_bound
    if a < 0:
        a = 0
    gamma_space = np.linspace(a, b, param_partition)
    #print "Gamma_space = ", gamma_space
    del(a, b)
    #-----------------

    print ("")
    print (" ----------------------------------------------------------------------------")
    print (" beta_space = "), beta_space
    print (" sigma_space = "), sigma_space
    print (" gamma_space = "), gamma_space
    print (" Number of points od the space for the current combined parameter space ="),\
        len(beta_space)*len(sigma_space)*len(gamma_space)
    print (" ----------------------------------------------------------------------------")
    print ("")
    ts_list = []
    ts_list_reduced = [] # samples sol. vector into corresponding obvserved ts
    param_distance = dict()  # dict to store param combination and distance
    R_squared = 0 # Coeffiecient of Determination

    for betas in beta_space:
        for sigmas in sigma_space:
            for gammas in gamma_space:
                print ("     beta, sigma, gamma: "), betas, sigmas, gammas
                my_sol = f_aggregated_model_solution(m_factor, H, Y0, t_start,
                                                     t_end, t_inc, t_range,
                                                     epi_classes,
                                                     betas, sigmas, gammas)
                my_sol_ts = f_aggregated_model_ts(my_sol)
                ts_list.append(my_sol_ts)

                # I sample the points in the solution that match in time with
                # the observed data so we can comapre them:
                my_sol_ts_reduced = f_sample_solution(my_observed_chik_ts,
                                                      my_sol_ts)
                ts_list_reduced.append(my_sol_ts_reduced)

                ## Here I calculate the distance between my_observed_chik_ts
                ## and my_sol_ts_reduced.
                distance = f_vector_distance(my_observed_chik_ts,
                                             my_sol_ts_reduced)

                # compute the coefficient of determination:
                R_squared = r_squared(my_observed_chik_ts, my_sol_ts_reduced)

                print ("")
                print ("")
                print ("     =========")
                print ("     distance: ", distance)
                print ("     R-squared: ", R_squared)

                # dictionary with keys=betas,sigmas,gammas,
                # and values=distance,R-squared:
                param_distance[(betas, sigmas, gammas)] = distance
                #print "Param_distance", param_distance
                print ("     =========")
                print ("")
                print ("")
    #for i in ts_list:
    #    i.plot(color='Crimson', linewidth=3, alpha=.45)
    #del(i)
    #pl.show()

    #my_observed_chik_ts.plot(color='DarkBlue', linewidth=4)
    #for i in ts_list_reduced:
    #    i.plot(color='Black', linewidth=3, alpha=.45)
    #del(i)
    #pl.show()

    # code to get the parameter set that gives the shortest distance between
    # observed data and model simulated data:
    min_distance = np.asarray(param_distance.values()).min() # shortest distance
    my_best_params = param_distance.keys()[param_distance.values().index(min_distance)]
    beta = my_best_params[0]
    sigma = my_best_params[1]
    gamma = my_best_params[2]
    print ('    beta = %f , sigma = %d, gamma = %d' %(beta, sigma, gamma))

    # plotting figs:
    # some configuration first:
    xlabelsize = 15
    ylabelsize = 15
    pl.rcParams['xtick.labelsize'] = xlabelsize
    pl.rcParams['ytick.labelsize'] = ylabelsize
    hfont = {'fontname':'Helvetica'}
    sefont = {'fontname':'Serif'}

    #pl.figure(figsize=(15,7))
    #my_observed_chik_ts.plot(color='Crimson', linewidth=5, alpha=.7)
    #my_last_plot = len(beta_space)*len(sigma_space)*len(gamma_space)
    #ts_list_reduced[my_last_plot - 1].plot(color = 'DarkGreen',
    #                                       linewidth=7, alpha=.7)
    #pl.title('Model fitting -least squares- Recursion #=%i, beta=%f, sigma=%f gamma=%f, distance between obs. and sol.=%f, init_bound=%f' %(recur, beta, sigma, gamma, min_distance, init_bound))
    #del(my_last_plot) # deleting temp variables
    #pl.ylabel("Number of infecteds", size=15, **sefont)
    #pl.xlabel("Time",size=15, **sefont)
    #pl.savefig('./figures/AGGREGATED_fitted_model_experiment=%f_recur=%d_beta=%f_sigma=%f_gamma=%f_partition=%d_distance=%f_initbound=%f.png'
    #           %(experiment, recur, beta, sigma, gamma, param_partition, min_distance, init_bound), dpi=400)

    ##pl.savefig('./figures/fitted_model_recur=%d_beta=%f_sigma=%f_gamma=%f_partition=%d_.jpg'
    #           %(recur, beta, sigma, gamma, param_partition), dpi=400)
    ##pl.savefig('fitted_model_beta%d_sigma%d_gamma%d_recur%d_partition%d_.jpg'
    #           %beta %sigma %gamma %recur %param_partition, dpi=400)
    #pl.show()
    del(experiment)

    print ("")
    print (' ------------------------------------------------------')
    print (" shortest distance is: "), min_distance
    print (" best parameter set is beta, sigma, gamma = "), my_best_params
    print (' ------------------------------------------------------')

    #end = time.time()
    #print ("----------")
    #print ("execution time: ", end - start, " secs")
    #print ("converted to minutes: ", (end - start)/60, "minutes")
    #print ("----------")
    del(beta, sigma, gamma)
    #return(param_distance)
    return(my_best_params) # returning best parameters beta, sigma, gamma

#################################################################################
# MINIMUM_SQUARE_ERROR2(MY_OBSERVED_CHIK_TS, H0, Y0, T_START, T_END, T_INC,
# T_RANGE,EP_CLASSES, BETA, SIGMA, GAMMA, BETA_SPACE, SIGMA_SPACE, GAMMA_SPACE)
def minimum_square_error2(my_observed_chik_ts,
                         H, Y0, t_start, t_end, t_inc, t_range, epi_classes,
                         init_beta, init_sigma, init_gamma,
                         init_bound, param_partition, recur, experiment,
                         beta_space, sigma_space, gamma_space):
    """
    This functions takes:
        - my_observed_chik_ts, a pandas time series object representing the
        aggregated beserved data of chik incidencein the Caribea outbreak
        reported by Tycho,
        - H, a networkx object that represents the oberved fliht connection for
        the area on those yeatas (2013) according openflight.org,
        - Y0, t_start, t_end, t_inc, t_range, conditions needed for teh ODE
        solver,
        - beta, sigma, and gamma , teh epi parameters to optimze by minimum mean
        squares
    This function gives back:
        - a tuple consisting of the  parameters (beta, sigma, gamma) that
        produces the solution with the shortet distance to the observed data
        for the corresponding discretization  of the theoretical solution
    """

    #start = time.time()

    print ("")
    print (" ----------------------------------------------------------------------------")
    print (" beta_space = "), beta_space
    print (" sigma_space = "), sigma_space
    print (" gamma_space = "), gamma_space
    print (" Number of points od the space for the current combined parameter space ="),\
        len(beta_space)*len(sigma_space)*len(gamma_space)
    print (" ----------------------------------------------------------------------------")
    print ("")
    ts_list = []
    ts_list_reduced = [] # samples sol. vector into corresponding obvserved ts
    param_distance = dict()  # dict to store param combination and distance
    R_squared = 0 # Coeffiecient of Determination

    for betas in beta_space:
        for sigmas in sigma_space:
            for gammas in gamma_space:
                print ("     beta, sigma, gamma: "), betas, sigmas, gammas
                my_sol = f_aggregated_model_solution(m_factor, H, Y0, t_start,
                                                     t_end, t_inc, t_range,
                                                     epi_classes,
                                                     betas, sigmas, gammas)
                my_sol_ts = f_aggregated_model_ts(my_sol)
                ts_list.append(my_sol_ts)

                # I sample the points in the solution that match in time with
                # the observed data so we can comapre them:
                my_sol_ts_reduced = f_sample_solution(my_observed_chik_ts,
                                                      my_sol_ts)
                ts_list_reduced.append(my_sol_ts_reduced)

                ## Here I calculate the distance between my_observed_chik_ts
                ## and my_sol_ts_reduced.
                distance = f_vector_distance(my_observed_chik_ts,
                                             my_sol_ts_reduced)

                # compute the coefficient of determination:
                R_squared = r_squared(my_observed_chik_ts, my_sol_ts_reduced)

                print ("")
                print ("")
                print ("     =========")
                print ("     distance: ", distance)
                print ("     R-squared: ", R_squared)

                # dictionary with keys=betas,sigmas,gammas,
                # and values=distance,R-squared:
                param_distance[(betas, sigmas, gammas)] = distance
                #print "Param_distance", param_distance
                print ("     =========")
                print ("")
                print ("")
    #for i in ts_list:
    #    i.plot(color='Crimson', linewidth=3, alpha=.45)
    #del(i)
    #pl.show()

    #my_observed_chik_ts.plot(color='DarkBlue', linewidth=4)
    #for i in ts_list_reduced:
    #    i.plot(color='Black', linewidth=3, alpha=.45)
    #del(i)
    #pl.show()

    # code to get the parameter set that gives the shortest distance between
    # observed data and model simulated data:
    min_distance = np.asarray(param_distance.values()).min() # shortest distance
    my_best_params = param_distance.keys()[param_distance.values().index(min_distance)]
    beta = my_best_params[0]
    sigma = my_best_params[1]
    gamma = my_best_params[2]
    print ('    beta = %f , sigma = %d, gamma = %d' %(beta, sigma, gamma))

    # plotting figs:
    # some configuration first:
    xlabelsize = 15
    ylabelsize = 15
    pl.rcParams['xtick.labelsize'] = xlabelsize
    pl.rcParams['ytick.labelsize'] = ylabelsize
    hfont = {'fontname':'Helvetica'}
    sefont = {'fontname':'Serif'}

    #pl.figure(figsize=(15,7))
    #my_observed_chik_ts.plot(color='Crimson', linewidth=5, alpha=.7)
    #my_last_plot = len(beta_space)*len(sigma_space)*len(gamma_space)
    #ts_list_reduced[my_last_plot - 1].plot(color = 'DarkGreen',
    #                                       linewidth=7, alpha=.7)
    #pl.title('Model fitting -least squares- Recursion #=%i, beta=%f, sigma=%f gamma=%f, distance between obs. and sol.=%f, init_bound=%f' %(recur, beta, sigma, gamma, min_distance, init_bound))
    #del(my_last_plot) # deleting temp variables
    #pl.ylabel("Number of infecteds", size=15, **sefont)
    #pl.xlabel("Time",size=15, **sefont)
    #pl.savefig('./figures/AGGREGATED_fitted_model_experiment=%f_recur=%d_beta=%f_sigma=%f_gamma=%f_partition=%d_distance=%f_initbound=%f.png'
    #           %(experiment, recur, beta, sigma, gamma, param_partition, min_distance, init_bound), dpi=400)

    ##pl.savefig('./figures/fitted_model_recur=%d_beta=%f_sigma=%f_gamma=%f_partition=%d_.jpg'
    #           %(recur, beta, sigma, gamma, param_partition), dpi=400)
    ##pl.savefig('fitted_model_beta%d_sigma%d_gamma%d_recur%d_partition%d_.jpg'
    #           %beta %sigma %gamma %recur %param_partition, dpi=400)
    #pl.show()
    del(experiment)

    print ("")
    print (' ------------------------------------------------------')
    print (" shortest distance is: "), min_distance
    print (" best parameter set is beta, sigma, gamma = "), my_best_params
    print (' ------------------------------------------------------')

    #end = time.time()
    #print ("----------")
    #print ("execution time: ", end - start, " secs")
    #print ("converted to minutes: ", (end - start)/60, "minutes")
    #print ("----------")
    del(beta, sigma, gamma)
    #return(param_distance)
    return(my_best_params) # returning best parameters beta, sigma, gamma


#################################################################################
# MINIMUM_SQUARE_ERROR_NODE(MY_OBSERVED_CHIK_TS, Y0, T_START, T_END, T_INC,
# T_RANGE,EP_CLASSES, BETA, SIGMA, GAMMA)
def minimum_square_error_node(my_observed_chik_ts, H,
                         Y0, t_start, t_end, t_inc, t_range, epi_classes,
                         init_beta, init_sigma, init_gamma,
                         init_bound, numberof_recurrences, param_partition, recur, experiment):
    """
    This functions takes:
        - my_observed_chik_ts, a pandas time series object representing the
        aggregated beserved data of chik incidencein the Caribea outbreak
        reported by Tycho,
        - Y0, t_start, t_end, t_inc, t_range, conditions needed for teh ODE
        solver,
        - beta, sigma, and gamma , teh epi parameters to optimze by minimum mean
        squares
    This function gives back:
        - a tuple consisting of the  parameters (beta, sigma, gamma) that
        produces the solution with the shortest distance to the observed data
        for the corresponding discretization  of the theoretical solution
    """

    #start = time.time()
    print("PARAMETER PARTITION = ", param_partition)

    a = init_beta - init_bound
    b = init_beta + init_bound
    if a < 0:
        a = 0
    beta_space = np.linspace(a, b, param_partition)
    #print "Beta_space = ", beta_space
    del (a,b)
    #-----------------

    a = init_sigma - init_bound
    b = init_sigma + init_bound
    if a < 0:
        a = 0
    sigma_space = np.linspace(a, b, param_partition)
    #print "Sigma_space = ", sigma_space
    del(a,b)
    #-----------------

    a = init_gamma - init_bound
    b = init_gamma + init_bound
    if a < 0:
        a = 0
    gamma_space = np.linspace(a, b, param_partition)
    #print "Gamma_space = ", gamma_space
    del(a, b)
    #-----------------

    print ("")
    print (" ----------------------------------------------------------------------------")
    print (" beta_space = "), beta_space
    print (" sigma_space = "), sigma_space
    print (" gamma_space = "), gamma_space
    print (" Number of points od the space for the current combined parameter space ="),\
        len(beta_space)*len(sigma_space)*len(gamma_space)
    print (" ----------------------------------------------------------------------------")
    print ("")
    ts_list = []
    ts_list_reduced = [] # samples sol. vector into corresponding obvserved ts
    param_distance = dict()  # dict to store param combination and distance
    R_squared = 0 # Coeffiecient of Determination

    for betas in beta_space:
        for sigmas in sigma_space:
            for gammas in gamma_space:
                print ("     beta, sigma, gamma: "), betas, sigmas, gammas
                my_sol = f_aggregated_model_solution(m_factor, H, Y0, t_start,
                                                     t_end, t_inc, t_range,
                                                     epi_classes,
                                                     betas, sigmas, gammas)
                my_sol_ts = f_aggregated_model_ts(my_sol)
                ts_list.append(my_sol_ts)

                # I sample the points in the solution that match in time with
                # the observed data so we can comapre them:
                my_sol_ts_reduced = f_sample_solution(my_observed_chik_ts,
                                                      my_sol_ts)
                ts_list_reduced.append(my_sol_ts_reduced)

                ## Here I calculate the distance between my_observed_chik_ts
                ## and my_sol_ts_reduced.
                distance = f_vector_distance(my_observed_chik_ts,
                                             my_sol_ts_reduced)

                # compute the coefficient of determination:
                R_squared = r_squared(my_observed_chik_ts, my_sol_ts_reduced)

                print ("")
                print ("")
                print ("     =========")
                print ("     distance: ", distance)
                print ("     R-squared: ", R_squared)

                # dictionary with keys=betas,sigmas,gammas,
                # and values=distance,R-squared:
                param_distance[(betas, sigmas, gammas)] = distance
                #print "Param_distance", param_distance
                print ("     =========")
                print ("")
                print ("")
    #for i in ts_list:
    #    i.plot(color='Crimson', linewidth=3, alpha=.45)
    #del(i)
    #pl.show()

    #my_observed_chik_ts.plot(color='DarkBlue', linewidth=4)
    #for i in ts_list_reduced:
    #    i.plot(color='Black', linewidth=3, alpha=.45)
    #del(i)
    #pl.show()

    # code to get the parameter set that gives the shortest distance between
    # observed data and model simulated data:
    min_distance = np.asarray(param_distance.values()).min() # shortest distance
    my_best_params = param_distance.keys()[param_distance.values().index(min_distance)]
    beta = my_best_params[0]
    sigma = my_best_params[1]
    gamma = my_best_params[2]
    print ('    beta = %f , sigma = %d, gamma = %d' %(beta, sigma, gamma))

    # plotting figs:
    # some configuration first:
    xlabelsize = 15
    ylabelsize = 15
    pl.rcParams['xtick.labelsize'] = xlabelsize
    pl.rcParams['ytick.labelsize'] = ylabelsize
    hfont = {'fontname':'Helvetica'}
    sefont = {'fontname':'Serif'}

    pl.figure(figsize=(15,7))
    my_observed_chik_ts.plot(color='Crimson', linewidth=5, alpha=.7)
    my_last_plot = len(beta_space)*len(sigma_space)*len(gamma_space)
    ts_list_reduced[my_last_plot - 1].plot(color = 'DarkGreen',
                                           linewidth=7, alpha=.7)
    pl.title('Model fitting -least squares- Recursion #=%i, beta=%f, sigma=%f gamma=%f, distance between obs. and sol.=%f, init_bound=%f' %(recur, beta, sigma, gamma, min_distance, init_bound))
    del(my_last_plot) # deleting temp variables
    pl.ylabel("Number of infecteds", size=15, **sefont)
    pl.xlabel("Time",size=15, **sefont)
    pl.savefig('./figures/fitted_model_experiment=%f_recur=%d_beta=%f_sigma=%f_gamma=%f_partition=%d_distance=%f_initbound=%f.png'
               %(experiment, recur, beta, sigma, gamma, param_partition, min_distance, init_bound), dpi=400)

    ##pl.savefig('./figures/fitted_model_recur=%d_beta=%f_sigma=%f_gamma=%f_partition=%d_.jpg'
    #           %(recur, beta, sigma, gamma, param_partition), dpi=400)
    ##pl.savefig('fitted_model_beta%d_sigma%d_gamma%d_recur%d_partition%d_.jpg'
    #           %beta %sigma %gamma %recur %param_partition, dpi=400)
    #pl.show()
    del(experiment)

    print ("")
    print (' ------------------------------------------------------')
    print (" shortest distance is: "), min_distance
    print (" best parameter set is beta, sigma, gamma = "), my_best_params
    print (' ------------------------------------------------------')

    #end = time.time()
    #print ("----------")
    #print ("execution time: ", end - start, " secs")
    #print ("converted to minutes: ", (end - start)/60, "minutes")
    #print ("----------")
    del(beta, sigma, gamma)
    #return(param_distance)
    return(my_best_params) # returning best parameters beta, sigma, gamma


#################################################################################
## FITTED_MODEL_NODE
def fitted_model_node(my_ts, H, Y0, t_start, t_end, t_inc, t_range, epi_classes,
                 init_beta, init_sigma, init_gamma, init_bound,
                 numberof_recurrences, param_partition, experiment):
    """
    this function receives:
        - Y0 (numpy vector), t_start, t_end, t_inc, t_range, epi_classes,
        - init_beta. init_sigma, init_gamma: initial values for the parameters
        to be fitted,
        - init_param_range: (tuple) initial range for parameters: usually (0,
        init_param_value + b)
        - numberof_recurrences: (int) this number account for the number of
        recurrences to fit the parameters: (int) number of recurrences over the
        parameter space to fit the parameters
        - param_partition: (int) partition used to explore recurrently the
        parameter space during the fitting process.
    this function gives back:
        - the best set of parameters, beta, sigma, gamma, which (discretized)
        simulation curves has the shortest distance to the observed chik data.
        This would be the least squared best fitted parameters.
    """
    print ('')
    print ('')
    print ('Total Number of Recurrences: ', numberof_recurrences)
    print ('=================================')
    print ('Fitting observed data to model in node: ', node)
    print ('==================================')
    print ('')
    print ('')

    # for fitting the initial stage of the epidemic let the following code,
    # these code lines shorten the initial obseved time series from the initial
    # day up to the 2014-05-11:
    my_ts.index = pd.DatetimeIndex(my_ts.index)
    my_ts = my_ts['2013-10-20':'2014-05-11']

    experiment = np.random.random_sample() # unique  naming the experiment ID

    for recur in np.arange(0, numberof_recurrences):
        print ('')
        print ('*')
        print ('*')
        print ('*')
        print ('-----------------------------')
        print ('-----------------------------')
        print ('-----------------------------')
        print ('------------------------------------------------------')
        print ("Recurrence #: "), recur
        print ('')
        my_fitted_params = minimum_square_error_node(my_ts, H, Y0,
                                                t_start, t_end, t_inc, t_range,
                                                epi_classes,
                                                init_beta, init_sigma, init_gamma,
                                                init_bound, numberof_recurrences, param_partition, recur, experiment)

        init_beta = my_fitted_params[0]
        init_sigma = my_fitted_params[1]
        init_gamma = my_fitted_params[2]
        #init_bound = init_bound/(param_partition**(recur+1))
        init_bound = init_bound/(param_partition)
        #param_partition = param_partition**(recur+1)

        print ("best fitted beta in this recurrence = "), my_fitted_params[0]
        print ("best fitted sigma in this recurrence = "), my_fitted_params[1]
        print ("best fitted gamma in this recurrence = "), my_fitted_params[2]
        print ("new beta partition goes from: "),init_beta - init_bound
        print ("to: "), init_beta + init_bound
        print ('------------------------------------------------------')
        print ('-----------------------------')
        print ('-----------------------------')
        print ('-----------------------------')
        print ('*')
        print ('*')
        print ('*')
        print ('')
        print ('')

    return (my_fitted_params)


#################################################################################
## FITTED_MODEL
def fitted_model(my_ts, H, Y0, t_start, t_end, t_inc, t_range, epi_classes,
                 init_beta, init_sigma, init_gamma, init_bound,
                 numberof_recurrences, param_partition):
    """
    this function receives:
        - H (networkx object), Y0 (numpy vector), t_start, t_end, t_inc, t_range, epi_classes,
        - init_beta. init_sigma, init_gamma: initial values for the parameters
        to be fitted,
        - init_param_range: (tuple) initial range for parameters: usually (0,
        init_param_value + b)
        - numberof_recurrences: (int) this number account for the number of
        recurrences to fit the parameters: (int) number of recurrences over the
        parameter space to fit the parameters
        - param_partition: (int) partition used to explore recurrently the
        parameter space during the fitting process.
    this function gives back:
        - the best set of parameters, beta, sigma, gamma, which (discretized)
        simulation curves has the shortest distance to the observed chik data.
        This would be the least squared best fitted parameters.
    """
    print ('')
    print ('')
    print ('Total Number of Recurrences: ', numberof_recurrences)
    print ('=================================')
    print ('Fitting observed data to model....')
    print ('==================================')
    print ('')
    print ('')

    # for fitting the initial stage of the epidemic let the following code,
    # these code lines shorten the initial obseved time series from the initial
    # day up to the 2014-05-11:
    my_ts.index = pd.DatetimeIndex(my_ts.index)
    my_ts = my_ts['2013-10-20':'2014-05-11']

    experiment = np.random.random_sample() # solely naming the experiment ID

    for recur in np.arange(0, numberof_recurrences):
        print ('')
        print ('*')
        print ('*')
        print ('*')
        print ('-----------------------------')
        print ('-----------------------------')
        print ('-----------------------------')
        print ('------------------------------------------------------')
        print ("Recurrence #: "), recur
        print ('')
        my_fitted_params = minimum_square_error(my_ts,H, Y0,
                                                t_start, t_end, t_inc, t_range,
                                                epi_classes,
                                                init_beta, init_sigma, init_gamma,
                                                init_bound, param_partition, recur, experiment)

        init_beta = my_fitted_params[0]
        init_sigma = my_fitted_params[1]
        init_gamma = my_fitted_params[2]
        #init_bound = init_bound/(param_partition**(recur+1))
        init_bound = init_bound/(param_partition)
        #param_partition = param_partition**(recur+1)

        print ("best fitted beta in this recurrence = "), my_fitted_params[0]
        print ("best fitted sigma in this recurrence = "), my_fitted_params[1]
        print ("best fitted gamma in this recurrence = "), my_fitted_params[2]
        print ("new beta partition goes from: "),init_beta - init_bound
        print ("to: "), init_beta + init_bound
        print ('------------------------------------------------------')
        print ('-----------------------------')
        print ('-----------------------------')
        print ('-----------------------------')
        print ('*')
        print ('*')
        print ('*')
        print ('')
        print ('')

    return (my_fitted_params)


#################################################################################
## FITTED_MODEL2
def fitted_model2(my_ts, H, Y0, t_start, t_end, t_inc, t_range, epi_classes,
                 init_beta, init_sigma, init_gamma, init_bound,
                 numberof_recurrences, param_partition):
    """
    this function receives:
        - H (networkx object), Y0 (numpy vector), t_start, t_end, t_inc, t_range, epi_classes,
        - init_beta. init_sigma, init_gamma: initial values for the parameters
        to be fitted,
        - init_param_range: (tuple) initial range for parameters: usually (0,
        init_param_value + b)
        - numberof_recurrences: (int) this number account for the number of
        recurrences to fit the parameters: (int) number of recurrences over the
        parameter space to fit the parameters
        - param_partition: (int) partition used to explore recurrently the
        parameter space during the fitting process.
    this function gives back:
        - the best set of parameters, beta, sigma, gamma, which (discretized)
        simulation curves has the shortest distance to the observed chik data.
        This would be the least squared best fitted parameters.
    """
    print ('')
    print ('')
    print ('Total Number of Recurrences: ', numberof_recurrences)
    print ('=================================')
    print ('Fitting observed data to model....')
    print ('==================================')
    print ('')
    print ('')

    # for fitting the initial stage of the epidemic let the following code,
    # these code lines shorten the initial obseved time series from the initial
    # day up to the 2014-05-11:
    my_ts.index = pd.DatetimeIndex(my_ts.index)
    my_ts = my_ts['2013-10-20':'2014-05-11']

    experiment = np.random.random_sample() # solely naming the experiment ID

    for recur in np.arange(0, numberof_recurrences):
        print ('')
        print ('*')
        print ('*')
        print ('*')
        print ('-----------------------------')
        print ('-----------------------------')
        print ('-----------------------------')
        print ('------------------------------------------------------')
        print ("Recurrence #: "), recur
        print ('')
        my_fitted_params = minimum_square_error2(my_ts,H, Y0,
                                                t_start, t_end, t_inc, t_range,
                                                epi_classes,
                                                init_beta, init_sigma, init_gamma,
                                                init_bound, param_partition, recur, experiment,
                                                 beta_space, sigma_space, gamma_space)

        init_beta = my_fitted_params[0]
        init_sigma = my_fitted_params[1]
        init_gamma = my_fitted_params[2]
        #init_bound = init_bound/(param_partition**(recur+1))
        init_bound = init_bound/(param_partition)
        #param_partition = param_partition**(recur+1)

        print ("best fitted beta in this recurrence = "), my_fitted_params[0]
        print ("best fitted sigma in this recurrence = "), my_fitted_params[1]
        print ("best fitted gamma in this recurrence = "), my_fitted_params[2]
        print ("new beta partition goes from: "),init_beta - init_bound
        print ("to: "), init_beta + init_bound
        print ('------------------------------------------------------')
        print ('-----------------------------')
        print ('-----------------------------')
        print ('-----------------------------')
        print ('*')
        print ('*')
        print ('*')
        print ('')
        print ('')

    return (my_fitted_params)


#################################################################################
## DEGREE HISTOGRAM
def degree_histogram(H):
    """
    this function receives:
        - a networkx object
    this function gives back:
        - a plot of the degree histogram of the networkx object
    """

    degree_sequence=sorted(nx.degree(H).values(),reverse=True) # degree sequence
    #print "Degree sequence", degree_sequence

    no_grid = {'grid': False} # seaborn plots without no grid
    pl.rc('axes', **no_grid)
    xlabelsize = 15
    ylabelsize = 15
    pl.rcParams['xtick.labelsize'] = xlabelsize
    pl.rcParams['ytick.labelsize'] = ylabelsize
    dmax = max(degree_sequence)
    #pl.figure()
    # plot degree rank
    pl.loglog(degree_sequence,'b-',marker='o',
              color='DarkBlue', alpha=0.6,
              markersize=10, linewidth=3,
              markeredgecolor='DarkBlue', markeredgewidth=2)

    hfont = {'fontname':'Helvetica'}
    sefont = {'fontname':'Serif'}

    # not to be confused with a DEGREE DSITRIBUTION Plot
    pl.title("Degree rank plot", size=20, **hfont)
    #pl.title("Degree distribution plot", size=20, **hfont)
    pl.ylabel("(log) Degree", size=15, **sefont)
    pl.xlabel("(log) Rank",size=15, **sefont)

    # draw graph in inset
    pl.axes([0.45,0.45,0.45,0.45])
    Gcc=sorted(nx.connected_component_subgraphs(H), key = len, reverse=True)[0]
    pos=nx.spring_layout(Gcc)
    pl.axis('off')
    nx.draw_networkx_nodes(Gcc,pos,node_size=100, alpha=0.6,  node_color='Crimson')
    nx.draw_networkx_edges(Gcc,pos,alpha=0.6, edge_color='Crimson', width=2)
    ##pl.savefig("./figures/degree_rank_plot.png")


    # plot degree distribution
    no_grid = {'grid': True} # seaborn plots back to default
    pl.rc('axes', **no_grid)

    pl.figure()
    pl.ylim(0,6.75)
    pl.hist(degree_sequence, color="Crimson", alpha=1.0)
    pl.title("Degree distribution", size=20, **hfont)
    #pl.title("Degree distribution plot", size=20, **hfont)
    pl.ylabel("Frequency", size=15, **sefont)
    pl.xlabel("Degree",size=15, **sefont)

    ##pl.savefig("./figures/degree_distribution_plot.png")
    pl.show()


#################################################################################
# F_PLOT_MODEL(M_FACTOR, H, T_START, T_END, T_INC, T_RANGE, EPI_CLASSES, BETA,
# SIGMA, GAMMA)
def f_plot_model(m_factor, H, Y0, t_start, t_end, t_inc, \
                 t_range, epi_classes, beta, sigma, gamma):
    """
    this function receives:
        - a adjacency matrix  object M (numpy matrix)
        - a proportional  rescaling constant for the matrix connection m_factor
        - a network connection modle H (networkx object)
        - the number of epidemic clases (including the class for the total
        number of individuals, e.g. an SEIR model has 5 classes: S,E,I,R, and N
        since the model needs to keep track on thi as well)
    thif function gives back:
        - a plot of the models solution
    """
    M = nx.to_numpy_matrix(H) # making sure M corresponds to H
    M = M*m_factor # rescale matrix
    patches = M.shape[0]

    # solving the eqs.:
    YSOL = spi.odeint(f_diff_eqs, Y0, t_range, args=(epi_classes, M, beta, sigma, gamma))
    YSOL = YSOL*.15 # only 15% are observed or recorded
    YSOL = YSOL*.3 # only 15% are observed or recorded

    # PLotting the total sum of the infection (theory):
    Total_Inf = np.zeros(len(t_range))
    count = 2
    for i in range(0,30):
        Total_Inf = Total_Inf + YSOL[:,count]
        count = count + 5
    del  count

    pl.figure(figsize=(15,7))
    pl.plot(t_range, Total_Inf, linewidth=6, alpha=.45, color='DarkBlue')
    pl.title('Theoretical curve of infected in all islands, beta=%f, sigma=%f, gamma=%f' %(beta, sigma, gamma), fontsize=20)
    pl.xlabel('time (days)', fontsize=20)
    pl.legend(loc=1)
    #pl.savefig('./figures/fitted_model_beta=%f_sigma=%f_gamma=%f.png'
    #           %(beta, sigma, gamma), dpi=400)
    return()


#################################################################################
# F_NODE_MODEL_SOLUTION(M_FACTOR, H, Y_0, T_START, T_END, T_INC, T_RANGE,
# EPI_CLASSES, BETA, SIGMA, GAMMA)
def f_node_model_solution(m_factor, H, Y0, t_start, t_end, t_inc, \
                          t_range, epi_classes, beta, sigma, gamma):
    """
    this function receives:
        - a proportional  rescaling constant for the matrix connection m_factor
        - a network connection modle H (networkx object)
        - a Y0 vector with initial condictions
        - the number of epidemic clases (including the class for the total
        number of individuals, e.g. an SEIR model has 5 classes: S,E,I,R, and N
        since the model needs to keep track on thi as well)
    this function gives back:
        - a numpy matrix where each column is a solution of the SEIRN network
        model
    """
    M = nx.to_numpy_matrix(H) # making sure M corresponds to H
    M = M*m_factor # rescale matrix

    # solving the eqs.:
    YSOL = spi.odeint(f_diff_eqs, Y0, t_range, args=(epi_classes, M, beta, sigma, gamma))
    YSOL = YSOL*.15 # only 15% are observed or recorded
    YSOL = YSOL*.3 # only 15% are observed or recorded

    # selecting the solution for infecteds:
    #Total_Inf = np.zeros(len(t_range))
    node_inf_list = list()
    count = 2
    for i in range(0,30):
        #Total_Inf = Total_Inf + YSOL[:,count]
        node_inf_list.append(YSOL[:,count])
        count = count + 5
    del  count

    #pl.figure(figsize=(15,7))
    #pl.plot(t_range, Total_Inf, linewidth=6, alpha=.45, color='DarkBlue')
    #pl.title('Theoretical curve of infected in all islands', fontsize=20)
    #pl.xlabel('time (days)', fontsize=20)
    #pl.legend(loc=1)

    return (node_inf_list)


#################################################################################
# F_AGGREGATED_MODEL_SOLUTION(M_FACTOR, H, Y_0, T_START, T_END, T_INC, T_RANGE,
# EPI_CLASSES, BETA, SIGMA, GAMMA)
def f_aggregated_model_solution(m_factor, H, Y0, t_start, t_end, t_inc, \
                          t_range, epi_classes, beta, sigma, gamma):
    """
    this function receives:
        - a proportional  rescaling constant for the matrix connection m_factor
        - a network connection modle H (networkx object)
        - a Y0 vector with initial condictions
        - the number of epidemic clases (including the class for the total
        number of individuals, e.g. an SEIR model has 5 classes: S,E,I,R, and N
        since the model needs to keep track on thi as well)
    this function gives back:
        - a numpy vector with the sum of all node-solutions of the network seir
        model
    """
    M = nx.to_numpy_matrix(H) # making sure M corresponds to H
    M = M*m_factor # rescale matrix

    # solving the eqs.:
    YSOL = spi.odeint(f_diff_eqs, Y0, t_range, args=(epi_classes, M, beta, sigma, gamma))
    YSOL = YSOL*.15 # only 15% are observed or recorded
    # re-scaling factor
    # rescale_factor = 0.3
    YSOL = YSOL*rescale_factor # new adjustment

    # PLotting the total sum of the infection (theory):
    Total_Inf = np.zeros(len(t_range))
    count = 2
    for i in range(0,30):
        Total_Inf = Total_Inf + YSOL[:,count]
        count = count + 5
    del  count

    #pl.figure(figsize=(15,7))
    #pl.plot(t_range, Total_Inf, linewidth=6, alpha=.45, color='DarkBlue')
    #pl.title('Theoretical curve of infected in all islands', fontsize=20)
    #pl.xlabel('time (days)', fontsize=20)
    #pl.legend(loc=1)

    return (Total_Inf)


#################################################################################
# F_NODE_MODEL_TS
def f_node_model_ts(my_model_solution_nodes):
    """
    this function receives:
        - a list where each element is a np.array object with a solution for
        infected per node.
    this function gives back:
        - a list with pandas time series object for each solution or
        - a  dictionary with node: infected time series
    """
    # list:
    node_inf_list = list()
    for i in range(0, len(my_model_solution_nodes)):
        node_inf_list.append(f_aggregated_model_ts(my_model_solution_nodes[i]))

    # dictionary:
    d = {}
    count =  0
    for i in H.nodes():
        d[i] = f_aggregated_model_ts(my_model_solution_nodes[count])
        count = count + 1

    # return list:
    #return(node_inf_list)

    # return dictionary:
    return(d)



#################################################################################
# F_AGGREGATED_MODEL_TS(MY_MODEL_SOLUTION)
def f_aggregated_model_ts(my_model_solution):
    """
    this function receives:
        - a numpy vector solution
    this function gives back:
        - a  pandas ts (Series) object that represent the vector solution
        model
    """
    # convert model solution to a time series:
    rng = pd.date_range('20-10-2013', periods=len(my_model_solution), freq='D')
    my_model_ts =  pd.Series(my_model_solution, index=rng)
    return (my_model_ts)


#################################################################################
# F_MODEL_SOLUTION_PLOT(MY_MODEL_SOLUTION)
def f_model_solution_plot(my_model_solution):
    """
    this function receives:
        - a numpy vector solution
    this function gives back:
        - a  plot of a pandas ts (Series) object that represent the vector
        solution model
    """
    # convert model solution to a time series:
    rng = pd.date_range('20-10-2013', periods=len(my_model_solution), freq='H')
    my_model_ts =  pd.Series(my_model_solution, index=rng)
    pl.figure(figsize=(15,7))
    my_model_ts.plot(linewidth=6, alpha=.85, color='DarkGreen')

    pl.title('Theoretical curve of infected in all islands, beta=%f, sigma=%f, gamma=%f' %(beta, sigma, gamma), fontsize=20)
    pl.xlabel('time (weeks)', fontsize=20)

    ## for ploting the fitted short time series:
    ##pl.savefig('./figures/fitted_model_TSshort_beta=%f_sigma=%f_gamma=%f.png'
    #           %(beta, sigma, gamma), dpi=400, orientation='landscape')
    ##pl.savefig('./figures/fitted_model_TSshort.png',
    #           dpi=400, orientation='landscape') # for including in ms

    # for plotting the fitted long time series:
    #pl.savefig('./figures/fitted_model_TSlong_beta=%f_sigma=%f_gamma=%f.png'
    #            %(beta, sigma, gamma), dpi=400, orientation='landscape')
    #pl.savefig('./figures/fitted_model_TSlong.png',
    #           orientation='landscape', dpi=400)  # for including in ms
    return ()

#################################################################################
# MODEL_TS_SHORT(MY_MODEL)
def model_ts_short(my_model):
    """
    - This function consumes:
        * a pd.Series object tha corresponds to the modeled time series.
    - This function gives back:
        * a pd.Series object corresponding to the shortened modeled time series
        that corresponds to the dates between 2013-10-20 and 2014-05-11, that is
        203 days (6.8 months) of the initial outbreak.
    """
    return(my_model['2013-10-20': '2014-05-11'])


#################################################################################
# PLOT_OBS_AND_MODEL(MY_TS, MY_MODEL)
def plot_obs_and_model(my_ts, my_model):
    my_ts.plot(linewidth=3, color="Crimson", alpha=.6)
    my_model.plot(linewidth=3, color="blue", alpha=.6)


#################################################################################
### function obs_ts_short
def obs_ts_short(my_ts):
    """
    - This function consumes:
        * a pd.Series object tha corresponds to the observed time series.
    - This function gives back:
        * a pd.Series object corresponding to the observed shortened time series
        that corresponds to the dates between 2013-10-20 and 2014-05-11, that is
        203 days (6.8 months) of the initial outbreak.
    """
    # for fitting the initial stage of the epidemic let the following code,
    # these code lines shorten the initial obseved time series from the initial
    # day up to the 2014-05-11:2014-05-11
    my_ts.index = pd.DatetimeIndex(my_ts.index)
    my_ts = my_ts['2013-10-20':'2014-05-11']
    return(my_ts)


#################################################################################
## plot_obs_and_model(my_ts, my_model):
def plot_obs_and_model(my_ts, my_model):
    """
    - This function consumes:
        * wot pd.Series object that corresponds to both the shortened (203 days)
        observed and modeled time series (of chik cases).
    - This function gives back:
        * a matplotlib plot of the two series simultaneously.
        """
    # Here I produce the figure of the shortened time series and solution -both as
    # pd.Series objects- out of  previously produced both series individually
    # plotting figs:
    # some configuration first:
    xlabelsize = 15
    ylabelsize = 15
    pl.rcParams['xtick.labelsize'] = xlabelsize
    pl.rcParams['ytick.labelsize'] = ylabelsize
    hfont = {'fontname': 'Helvetica'}
    sefont = {'fontname': 'Serif'}
    pl.figure(figsize=(12, 6))
    my_obs_time_series = obs_ts_short(my_observed_chik_ts)
    my_obs_time_series.plot(linewidth=4, color="Crimson", alpha=.6)
    my_model_time_series = model_ts_short(my_model_ts)
    my_model_time_series.plot(linewidth=4, color="blue", alpha=.6)
    pl.title("Observed aggregated data of chikungunya cases in the region (red) and fitted SIR Network model (blue)", size=15)
    pl.xlabel("Time", size=15, **sefont)
    pl.ylabel("Number of infecteds", size=15, **sefont)



####################### Scenarios ###############################################
### Functions defined especifically for the requested scenarios #################
#################################################################################
#### Plot Scenario (A)
def f_plot_scenario_A(reps, H, ):
    """
    this function receives:
        - reps, number of experiments, repetitions
        - H, a networkx object with the connection flight model
    this function gives back:
        - a plot of the soultions of  scenario A converted into pd.Series
        objects indexed my time
    """
    time_series = list() # to store time series
    stats = list()   # to store stats (max number of infcss, time to max, etc)
    pl.figure(figsize=(15,7))

    # matrix to store solition:
    sol_matrix = np.zeros([t_range.shape[0], reps])

    for i in np.arange(0, reps):

        print("This is i = ", i)

        new_H = f_H_scenario_A(H)
        my_scenario_A_solution = f_aggregated_model_solution(m_factor, new_H,
                                                             Y0,t_start, t_end,
                                                             t_inc, t_range,
                                                             epi_classes, beta,
                                                             sigma, gamma)
        sol_matrix[:,i] = my_scenario_A_solution
        #my_scenario_A_plot = f_model_solution_plot(my_scenario_A_solution)
        my_scenario_A_ts = f_aggregated_model_ts(my_scenario_A_solution)
        #a list with my neighbours during this rep:
        my_neighbours = new_H.edge['Saint Martin'].keys()
        my_risky_pop = np.sum([my_new_H_scenario_A.node[neigh]['pop_size']
                         for neigh in my_neighbours])
        print("Closest Neighbour Risk Population, Scenario A: ", my_risky_pop)

        # a list for the pops in risk during this rep:
        time_series.append(my_scenario_A_ts)
        # store: #experiment, max of the series,
        stats.append([i, my_scenario_A_ts.max(),
                      t_range[my_scenario_A_solution.argmax()], my_risky_pop ])
        #del(new_H, my_scenario_A_solution, my_scenario_A_ts,
        #    my_neighbours, my_risky_pop)

    for i in time_series:
        i.plot(color='Crimson', linewidth=3, alpha=.45)
    del(i)

    # I plot the unmodded model for comparison:
    my_model_ts.plot(color='DarkBlue', linewidth=4)
    #pl.title('Scenario A: randomized connections for initial infected island, red curves.\n Blue curve is the model with the observed connection network ', fontsize=20)
    pl.title('Epidemic scenario 1: randomized connections for initial infected island, red curves.\n Blue curve is the model with the observed connection network ', fontsize=20)
    #pl.title('Scenario A: randomized connections for initial infected island, red curves.'
    #         + r' Blue curve is the model with the observed connection network ', fontsize=20)
    pl.xlabel('time (days)', fontsize=20)
    #pl.legend(loc=1)
    #pl.savefig("randomized_connections_saint_martin_scenario_A.png",dpi=400)
    #print(stats)
    #pl.show()

    # histograms of max peak and time to max peak of the experiment:
    histogram_maxatpeak = np.zeros(len(stats), dtype=float)
    histogram_timetomax = np.zeros(len(stats), dtype=float)
    histogram_popatrisk = np.zeros(len(stats), dtype=float)
    for i in range(0, len(stats)):
        histogram_maxatpeak[i] = stats[i][1]
        histogram_timetomax[i] = stats[i][2]
        histogram_popatrisk[i] = np.log(stats[i][3])
        if  stats[i][3] == 0:
            histogram_popatrisk[i] = 0
    del(i)

    pl.figure()
    pl.hist(histogram_maxatpeak, color='crimson', bins=10)
    pl.xlabel('Number of infected at maximum epidemic peak', fontsize=20)
    pl.ylabel('frequency', fontsize=20)
    #pl.savefig("infected_at_maximum_epidemic_peak_scenario_A.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist(histogram_timetomax, bins=10)
    pl.xlabel('Time in days to reach epidemic peak', fontsize=20)
    pl.ylabel('frequency', fontsize=20)
    #pl.savefig("time_to_peak_scenario_A.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist(histogram_popatrisk, bins=10)
    pl.xlabel('Population at risk on first degree', fontsize=20)
    pl.ylabel('frequency', fontsize=20)
    #pl.savefig("riksy_pop_scenario_A.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist2d(histogram_maxatpeak, histogram_timetomax, bins=15)
    pl.xlabel('Number of infected at maximum epidemic peak', fontsize=15)
    pl.ylabel('Time in days to reach epidmic peak', fontsize=15)
    pl.title('2D Histogram', fontsize=22)
    pl.colorbar()
    #pl.savefig("Hist1_scenario_A.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist2d(histogram_maxatpeak, histogram_popatrisk, bins=15)
    pl.xlabel('Number of infected at maximum epidemic peak', fontsize=15)
    pl.ylabel('Population at risk on first degree', fontsize=15)
    pl.title('2D Histogram', fontsize=22)
    pl.colorbar()
    #pl.savefig("Hist2_scenario_A.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist2d(histogram_timetomax, histogram_popatrisk, bins=15)
    pl.xlabel('Time in days to reach epidmic peak', fontsize=15)
    pl.ylabel('Population at risk on first degree', fontsize=15)
    pl.title('2D Histogram', fontsize=22)
    pl.colorbar()
    #pl.savefig("Hist3_scenario_A.png",dpi=400)
    #pl.show()

    # 3d-histogram (pending):
    return (sol_matrix)


#################################################################################
# F_SCENARIO_A(H)
def f_H_scenario_A(H):
    """
    this function receives:
        - a networkx object H, representing the flight connection model
    this function gives back:
        - A new modified copy of H acoording scenario (A)
    """
    #############################################################################
    # Scenario A)
    # Randomization of the connection of the first infected island -Saint Martin-
    # In this scenario I initially seed one infected in Saint Martin and for each
    # run I take Saint Martin's original links and randomly reconnect them to
    #############################################################################
    ########
    def f_new_H(H, weights):
        """
        this function receives:
            - the model network model ,
            - a list of link weights
        and gives back:
            - a copy of H but with randomized connections for Saint Martin according
            scenario A
        """
        my_random_H = H.copy()
        # I remove all the links off Saint Martin:
        for i in H.edge['Saint Martin'].keys():
            my_random_H.remove_edge('Saint Martin', i)
        del i

        mylist = H.nodes()
        mylist.remove('Saint Martin')

        # it randomly choses a weight out of the list
        my_random_weight = np.random.choice(weights, 1)[0]
        # it randomly choses a country out of the list
        np.random.choice(mylist, 1)[0]
        for i in weights:
            my_random_H.add_edge('Saint Martin', np.random.choice(mylist, 1)[0],
                                 weight = np.random.choice(weights, 1)[0])
        del mylist, my_random_weight

        return(my_random_H)
    ########

    # some smart code to get the edge weights out of the Saint Martin node in the
    # networkx object H:
    list_of_values = [H.edge[item].values()[0].values()
                      for item in H.edge['Saint Martin']]
    weights = [i[0] for i in list_of_values]
    new_H = f_new_H(H, weights) # new modded H for scenario (A)

    del (list_of_values, weights)
    return(new_H)


#################################################################################
# F_SHUFFLE_Y0_SCENARIO_B(H)
def f_shuffle_Y0_scenario_B(H):
    """
    this function receives:
        - a netowrkx object that represents the oberved connection flights
    and gives back:
        - a dictionary where the keys() contents the name of the initial
        infected node node and the values() contains the corresponding Y0
        for this particular initial condition needed for ode solver to
        simulate scenario B.
    Scneario B (scenarion 4 in the paper): Initial infected node chosen randomly
    and original network based on observed flightsis preserved
    Scenerio B is EPIDEMIC SCENARIO 4 in the paper
    Initial infected node randomly chosen, observed flight network preserved
    """
    patches = H.number_of_nodes()*epi_classes # #-of patches
    Y0 = np.zeros(patches)
    # fill the initial vector:
    N0 = np.zeros(shape=(H.number_of_nodes()), dtype=float)

    # fill N0:
    count = 0
    for i in H.nodes():
        # N0[count] = H.node[i]['density'] # for densities
        N0[count] = H.node[i]['pop_size'] # for pop_sizes
        count = count + 1
    del count

    # fill initial Y0 (without seeding it yet):
    count = 0
    for n in N0:
        Y0[count] = n # for class S
        Y0[count+4] = n # for class N
        count = count + 5
    del count, N0, n
    to_chose = [i for i in np.arange(0, 150, 5)]
    chosen  = np.random.choice(to_chose, 1)[0]

    # find the corresponding island country:
    for i in H.node:
        if H.node[i]['pop_size'] == Y0[chosen]:
            #print "HELL YEAH!"
            mychosen = i
            #print chosen, mychosen
            del(i)

    Y0[chosen] = Y0[chosen] - 1  # S = S - 1
    Y0[chosen + 2] = Y0[chosen + 2] + 1  # I = I + 1
    return({mychosen: Y0})


#################################################################################
# F_PLOT_SCENARIO_B(REPS, H, MY_RANDOMIZED_Y0_SCENARIO_B)
def f_plot_scenario_B(reps, H, my_randomized_Y0_scenario_B):
    """
    this function receives:
        - reps, number of experiments, repetitions
        - H, a networkx object with the connection flight model
        - my_randomized_Y0_scenario_B which is a dict with the key having the
        name of the initial infected node and the value is the corresponding
        initial Y0
    this function gives back:
        - a plot of the soultions of  scenario B converted into pd.Series
        objects indexed my time
    Scenerio B is EPIDEMIC SCENARIO 4 in the paper
    Initial infected node randomly chosen, observed flight network preserved
    """
    time_series = list() # to store time series
    stats = list()   # to store stats (max number of infcss, time to max, etc)
    pl.figure(figsize=(15,7))

    for i in np.arange(0, reps):  # reps
        # converting dict's values() to an one-dimension np array:
        randomized_Y0 = f_shuffle_Y0_scenario_B(H)
        # getting the initial infected node's name as a string:
        my_Y0 = np.asarray(randomized_Y0.values())[0]
        mynode = randomized_Y0.keys()[0]
        my_scenario_B_solution = f_aggregated_model_solution(m_factor, H,
                                                             my_Y0,t_start, t_end,
                                                             t_inc, t_range,
                                                             epi_classes, beta,
                                                             sigma, gamma)
        #my_scenario_B_plot = f_model_solution_plot(my_scenario_B_solution)
        my_scenario_B_ts = f_aggregated_model_ts(my_scenario_B_solution)
        #a list with my neighbours during this rep:
        my_neighbours = H.edge[mynode].keys()
        my_risky_pop = np.sum([H.node[neigh]['pop_size']
                               for neigh in my_neighbours])
        print("Starting infected node, Scenario B: ",mynode)
        print("Closest Neighbour Risk Population, Scenario B: ", my_risky_pop)
        print ("")
        # a list for the pops in risk during this rep:
        time_series.append(my_scenario_B_ts)
        # store: #experiment, max of the series,
        stats.append([i, my_scenario_B_ts.max(),
                      t_range[my_scenario_B_solution.argmax()], my_risky_pop ])
        #del(new_H, my_scenario_A_solution, my_scenario_A_ts,
        #    my_neighbours, my_risky_pop)
    del my_Y0, mynode, my_neighbours, my_risky_pop
    for i in time_series:
        i.plot(color='Crimson', linewidth=4, alpha=.45)
    del(i)

    # I plot the unmodded model for comparison:
    my_model_ts.plot(color='DarkBlue', linewidth=6)
    pl.title('Scenario B:Initial infected node randomly chosen, observed flight network preserved (red curves).\n Blue curve is the model with Saint Martin as initial infected node.', fontsize=20)
    #pl.title('Initial infected node randomly chosen, observed flight network preserved (red curves).'
    #         + r' Blue curve is the model with Saint Martin as initial infected
    #         node.', fontsize=20)
    pl.xlabel('time (days)', fontsize=20)
    #pl.legend(loc=1)
    #pl.savefig("randomized_initial_seed_scenario_B1.png",dpi=400)
    #print(stats)
    #pl.show()

    # histograms of max peak and time to max peak of the experiment:
    histogram_maxatpeak = np.zeros(len(stats), dtype=float)
    histogram_timetomax = np.zeros(len(stats), dtype=float)
    histogram_popatrisk = np.zeros(len(stats), dtype=float)
    for i in range(0, len(stats)):
        histogram_maxatpeak[i] = stats[i][1]
        histogram_timetomax[i] = stats[i][2]
        histogram_popatrisk[i] = np.log(stats[i][3])
        if  stats[i][3] == 0:
            histogram_popatrisk[i] = 0
    del(i)

    pl.figure()
    pl.hist(histogram_maxatpeak, color='crimson', bins=10)
    pl.xlabel('Number of infected at maximum epidemic peak', fontsize=20)
    pl.ylabel('frequency', fontsize=20)
    #pl.show()
    #pl.savefig("infected_at_maximum_epidemic_peak_scenario_B1.png",dpi=400)
    pl.figure()
    pl.hist(histogram_timetomax, bins=10)
    pl.xlabel('Time in days to reach epidemic peak', fontsize=20)
    pl.ylabel('frequency', fontsize=20)
    #pl.savefig("time_to_peak_scenario_B1.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist(histogram_popatrisk, bins=10)
    pl.xlabel('Population at risk on first degree', fontsize=20)
    pl.ylabel('frequency', fontsize=20)
    #pl.savefig("riksy_pop_scenario_B1.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist2d(histogram_maxatpeak, histogram_timetomax, bins=15)
    pl.xlabel('Number of infected at maximum epidemic peak', fontsize=15)
    pl.ylabel('Time in days to reach epidmic peak', fontsize=15)
    pl.title('2D Histogram - Scenario B', fontsize=22)
    pl.colorbar()
    #pl.savefig("Hist1_scenario_B1.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist2d(histogram_maxatpeak, histogram_popatrisk, bins=15)
    pl.xlabel('Number of infected at maximum epidemic peak', fontsize=15)
    pl.ylabel('Population at risk on first degree', fontsize=15)
    pl.title('2D Histogram - Scenario B', fontsize=22)
    pl.colorbar()
    #pl.savefig("Hist2_scenario_B1.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist2d(histogram_timetomax, histogram_popatrisk, bins=15)
    pl.xlabel('Time in days to reach epidmic peak', fontsize=15)
    pl.ylabel('Population at risk on first degree', fontsize=15)
    pl.title('2D Histogram - Scenario B', fontsize=22)
    #pl.savefig("Hist3_scenario_B1.png",dpi=400)
    pl.colorbar()
    #pl.show()

    # 3d-histogram (pending):
    return ()


#################################################################################
# F_PLOT_SCENARIO_B2(REPS, H, MY_RANDOMIZED_Y0_SCENARIO_B)
def f_plot_scenario_B2(H):
    """
    This is version 2  of f_plot_scenario_B
    this function receives:
        - H, a networkx object with the connection flight model
        name of the initial infected node and the value is the corresponding
        initial Y0
    this function gives back:
        - a plot of the solutions of  scenario B converted into pd.Series
        objects indexed my time
    This is Scenario 4 in the paper
    """
    time_series = list() # to store time series
    stats = list()   # to store stats (max number of infcss, time to max, etc)
    pl.figure(figsize=(15,7))
    Y0_position = 0

    # matrix to store solution:
    sol_matrix = np.zeros([t_range.shape[0], H.number_of_nodes()])
    num_exp = 0 # counter for experiments

    for mynode in H.nodes():
        #my_Y0 = np.asarray(randomized_Y0.values())[0]
        #I need my Y0 corresponding to mynode:

        patches = H.number_of_nodes()*epi_classes # #-of patches
        Y0 = np.zeros(patches)
        N0 = np.zeros(shape=(H.number_of_nodes()), dtype=float)

        count = 0
        for i in H.nodes():
            # N0[count] = H.node[i]['density'] # for densities
            N0[count] = H.node[i]['pop_size'] # for pop_sizes
            count = count + 1
        del count

        # filling the initial vector (without seeding it yet):
        count = 0
        for n in N0:
            Y0[count] = n # for class S
            Y0[count+4] = n # for class N
            count = count + 5
        del count, N0, n

        Y0[Y0_position + 2] = 1
        Y0[Y0_position + 4] = Y0[Y0_position + 4] - 1
        Y0_position = Y0_position + 5

        my_Y0 = Y0 # to be consequent with notation
        del(Y0)

        my_scenario_B_solution = f_aggregated_model_solution(m_factor, H,
                                                             my_Y0,t_start, t_end,
                                                             t_inc, t_range,
                                                             epi_classes, beta,
                                                             sigma, gamma)

        sol_matrix[:,num_exp] = my_scenario_B_solution
        print("Node: ",mynode )
        print(my_scenario_B_solution.shape)
        #my_scenario_B_plot = f_model_solution_plot(my_scenario_B_solution)
        my_scenario_B_ts = f_aggregated_model_ts(my_scenario_B_solution)
        #a list with my neighbours during this rep:
        my_neighbours = H.edge[mynode].keys()
        my_risky_pop = np.sum([H.node[neigh]['pop_size']
                               for neigh in my_neighbours])
        print("Starting infected node, Scenario B: ",mynode)
        print("Closest Neighbour Risk Population, Scenario B: ", my_risky_pop)
        print ("")
        # a list for the pops in risk during this rep:
        time_series.append(my_scenario_B_ts)
        # store: #experiment, max of the series,
        stats.append([i, my_scenario_B_ts.max(),
                      t_range[my_scenario_B_solution.argmax()], my_risky_pop ])
        #del(new_H, my_scenario_A_solution, my_scenario_A_ts,
        #    my_neighbours, my_risky_pop)
        num_exp = num_exp + 1
    del my_Y0, mynode, my_neighbours, my_risky_pop
    for i in time_series:
        i.plot(color='Crimson', linewidth=4, alpha=.65)
    del(i)

    # I plot the unmodded model for comparison:
    my_model_ts.plot(color='DarkBlue', linewidth=6)
    #pl.title('Scenario B: Initial infected node randomly chosen, observed flight network is preserved (red curves).\n Blue curve is the model with Saint Martin as intial infected node.', fontsize=20)
    pl.title('Epidemic scenario four: Initial infected node randomly chosen, observed flight network is preserved (red curves).\n Blue curve is the model with Saint Martin as intial infected node.', fontsize=20)
    #pl.title('Initial infected node randomly chosen, observed flight network is preserved (red curves).'
    #         + r' Blue curve is the model with Saint Martin as intial infected node.', fontsize=20)
    pl.xlabel('time (days)', fontsize=20)
    #pl.legend(loc=1)
    #pl.savefig("randomized_initial_seed_scenario_B.png",dpi=400)
    #pl.savefig("randomized_initial_seed_scenario_4.png",dpi=400)
    #print(stats)
    #pl.show()

    # histograms of max peak and time to max peak of the experiment:
    histogram_maxatpeak = np.zeros(len(stats), dtype=float)
    histogram_timetomax = np.zeros(len(stats), dtype=float)
    histogram_popatrisk = np.zeros(len(stats), dtype=float)
    for i in range(0, len(stats)):
        histogram_maxatpeak[i] = stats[i][1]
        histogram_timetomax[i] = stats[i][2]
        histogram_popatrisk[i] = np.log(stats[i][3])
        if stats[i][3] == 0:
            histogram_popatrisk[i] = 0 # to correct the problem with -Inf

    del(i)

    pl.figure()
    pl.hist(histogram_maxatpeak, color='crimson', bins=10)
    pl.xlabel('Number of infected at maximum epidemic peak', fontsize=20)
    pl.ylabel('frequency', fontsize=20)
    #pl.show()
    #pl.savefig("infected_at_maximum_epidemic_peak_scenario_B.png",dpi=400)
    pl.figure()
    pl.hist(histogram_timetomax, bins=10)
    pl.xlabel('Time in days to reach epidemic peak', fontsize=20)
    pl.ylabel('frequency', fontsize=20)
    #pl.savefig("time_to_peak_scenario_B.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist(histogram_popatrisk, bins=10)
    pl.xlabel('Population at risk on first degree', fontsize=20)
    pl.ylabel('frequency', fontsize=20)
    #pl.savefig("riksy_pop_scenario_B.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist2d(histogram_maxatpeak, histogram_timetomax, bins=15)
    pl.xlabel('Number of infected at maximum epidemic peak', fontsize=15)
    pl.ylabel('Time in days to reach epidmic peak', fontsize=15)
    pl.title('2D Histogram - Scenario B', fontsize=22)
    pl.colorbar()
    #pl.savefig("Hist1_scenario_B.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist2d(histogram_maxatpeak, histogram_popatrisk, bins=15)
    pl.xlabel('Number of infected at maximum epidemic peak', fontsize=15)
    pl.ylabel('Population at risk on first degree', fontsize=15)
    pl.title('2D Histogram - Scenario B', fontsize=22)
    pl.colorbar()
    #pl.savefig("Hist2_scenario_B.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist2d(histogram_timetomax, histogram_popatrisk, bins=15)
    pl.xlabel('Time in days to reach epidmic peak', fontsize=15)
    pl.ylabel('Population at risk on first degree', fontsize=15)
    pl.title('2D Histogram - Scenario B', fontsize=22)
    #pl.savefig("Hist3_scenario_B.png",dpi=400)
    pl.colorbar()
    #pl.show()

    # 3d-histogram (pending):
    return (sol_matrix)


#################################################################################
# F_H_SCENARIO_C(H)
def f_H_scenario_C(H):
    """
    this function receives:
        - a networkx object H, representing the flihgt connection model
    this function gives back:
        - A new modified copy of H acoording scenario (C)
    """
    #############################################################################
    # Scenario C)
    # Connect first infected island (Saint Martin) to all others and make a run
    ########
    my_connected_H = H.copy()
    # I remove all the links off Saint Martin:
    for i in H.edge['Saint Martin'].keys():
        my_connected_H.remove_edge('Saint Martin', i)
    del i

    #############################################################################
    def f_new_H(H, weights):
        """
        this function receives:
            - the model network model ,
            - a list of link weights
        and gives back:
            - a copy of H but with Saint Martin connected to all others nodes
            according scenario C
        """
        my_connected_H = H.copy()
        # I remove all the links off Saint Martin:
        for i in H.edge['Saint Martin'].keys():
            my_connected_H.remove_edge('Saint Martin', i)
        del i

        mylist = H.nodes()
        mylist.remove('Saint Martin')

        # loop over all other nodes:
        # it randomly choses a weight out of the list
        for i in mylist:
            my_connected_H.add_edge('Saint Martin', i, weight =
                                    np.random.choice(weights, 1)[0])
        del mylist

        return(my_connected_H)
    ########

    # some smart code to get the edge weights out of the Saint Martin node in the
    # networkx object H:
    list_of_values = [H.edge[item].values()[0].values()
                      for item in H.edge['Saint Martin']]
    weights = [i[0] for i in list_of_values]
    new_H = f_new_H(H, weights) # new modded H for scenario (A)

    del (list_of_values, weights)
    return(new_H)


#################################################################################
#### Plot Scenario (C)
def f_plot_scenario_C(H, Y0):
    """
    this function receives:
        - H, a networkx object with the connection flight model according
        scenario C
    theis function gives back:
        - a plot of the soultions of  scenario C converted into pd.Series
        objects indexed my time
    """
    time_series = list() # to store time series
    stats = list()   # to store stats (max number of infcss, time to max, etc)
    pl.figure(figsize=(15,7))

    new_H = f_H_scenario_C(H)
    my_scenario_C_solution = f_aggregated_model_solution(m_factor, new_H,
                                                             Y0,t_start, t_end,
                                                             t_inc, t_range,
                                                             epi_classes, beta,
                                                             sigma, gamma)
    #my_scenario_C_plot = f_model_solution_plot(my_scenario_C_solution)
    my_scenario_C_ts = f_aggregated_model_ts(my_scenario_C_solution)
    #a list with my neighbours during this rep:
    my_neighbours = new_H.edge['Saint Martin'].keys()
    my_risky_pop = np.sum([my_new_H_scenario_C.node[neigh]['pop_size']
                           for neigh in my_neighbours])
    print("Closest Neighbour Risk Population, Scenario C: ", my_risky_pop)

    # a list for the pops in risk during this rep:
    time_series.append(my_scenario_C_ts)
    # store:  max of the series,
    stats.append([ my_scenario_C_ts.max(),
                  t_range[my_scenario_C_solution.argmax()], my_risky_pop ])
    #del(new_H, my_scenario_A_solution, my_scenario_A_ts,
    #    my_neighbours, my_risky_pop)

    # plotting figs:
    # some configuration first:
    xlabelsize = 12
    ylabelsize = 12
    pl.rcParams['xtick.labelsize'] = xlabelsize
    pl.rcParams['ytick.labelsize'] = ylabelsize
    hfont = {'fontname':'Helvetica'}
    sefont = {'fontname':'Serif'}
    my_scenario_C_ts.plot(color='Crimson', linewidth=6, alpha=.85)

    # I plot the unmodded model for comparison:
    my_model_ts.plot(color='DarkBlue', linewidth=8)
    #pl.title('Scenario C: Fully connected initial infected island, red curves.\n Blue curve is the model with the observed connection network ', fontsize=20)
    pl.title('Epidemic scenario three: Fully connected initial infected island, red curves.\n Blue curve is the model with the observed connection network ', fontsize=20,)# **sefont)
    pl.xlabel('time (days)', fontsize=20,)# **sefont)
    pl.ylabel('cases', fontsize=20,)# **sefont)
    #pl.legend(loc=1)
    #pl.savefig("fully_connected_saint_martin_scenario_C.png",dpi=400)
    #print(stats)
    #pl.show()

    # 3d-histogram (pending):
    return (my_scenario_C_ts)


#################################################################################
# F_H_SCENARIO_D(H, reps)
def f_H_scenario_D(H):
    """
    this function receives a networkx object representing a network of observed
    flight connections and it gives back an copy of the netwok based on the
    scenario D.
    """
    ## Randomize the network keeping the same number of links and
    ## nodes:
    # I use an uncorrelated random Erdos-Renyi G(n,p) graph:
    # where n is the number of nodes and
    # p is the probability of link formation
    # so n nodes are connected through l edges wich are chosen
    # randomly from n(n-1)/2 possible configurations
    # Every pair of nodes are connected with probability p
    # The total number of edges is a random variable
    # with expected value p*n(n-1)/2
    # hence p = 2l/(n(n-1)).
    # Define then a function for generating a random network
    # with links H.number_of_edges() and with
    # nodes H.number_of_nodes. This would create random networks models with probabilities
    # corresponding to the observed nodes and number of links.
    p = 2.*H.number_of_edges() / (H.number_of_nodes()*(H.number_of_nodes()) - 1.)
    my_random_network = nx.erdos_renyi_graph(H.number_of_nodes(), p)
    # relabel nodes:
    key = my_random_network.nodes()
    value = H.nodes()
    key_and_value = zip(key, value)
    # my dict:
    mapping = {}
    for key, value in key_and_value:
        mapping[key] = value
    del (key, value, key_and_value)
    my_random_network = nx.relabel_nodes(my_random_network, mapping)
    # add the corresponding weights
    # so I need to recover weights from the oroginal network/adjacency matrix:
    M = nx.to_numpy_matrix(H)
    # set of existing weight values:
    weights = np.unique(np.asarray(M))

    # assign to each of the randomly formed link and random weight out of
    # existing weights
    for i in my_random_network.edges():
        my_random_network.add_edge(i[0],i[1],
                                   weight = np.random.choice(weights, 1))

    # assigning pop_sizes:
    for i in H.nodes():
        my_random_network.node[i]['pop_size'] = H.node[i]['pop_size']

    return(my_random_network)


#################################################################################
#### Plot Scenario (D)
def f_plot_scenario_D(reps, H, Y0):
    """
    this function receives:
        - reps, number of experiments, repetitions
        - H, a networkx object with the observed connection flight model
        - ramdom_H. a networkx object with connection flight model according to
        scenario D,
        - Y0, and array  for the initial condition for the ODE
    theis functio gives back:
        - a plot of the soultions of  scenario D  converted into pd.Series
        objects indexed my time
        - a (time x reps) matrix where columns are repetitions of the
        experimets' solutions
    """
    time_series = list() # to store time series
    stats = list()   # to store stats (max number of infcss, time to max, etc)
    pl.figure(figsize=(15,7))

    # matrix to store solition:
    sol_matrix = np.zeros([t_range.shape[0], reps])

    for i in np.arange(0, reps):
        new_H = f_H_scenario_D(H)
        my_scenario_D_solution = f_aggregated_model_solution(m_factor, new_H,
                                                             Y0,t_start, t_end,
                                                             t_inc, t_range,
                                                             epi_classes, beta,
                                                             sigma, gamma)
        sol_matrix[:,i] = my_scenario_D_solution
        #my_scenario_D_plot = f_model_solution_plot(my_scenario_D_solution)
        my_scenario_D_ts = f_aggregated_model_ts(my_scenario_D_solution)
        #a list with my neighbours during this rep:
        my_neighbours = new_H.edge['Saint Martin'].keys()
        my_risky_pop = np.sum([my_new_H_scenario_D.node[neigh]['pop_size']
                         for neigh in my_neighbours])
        print("Closest Neighbour Risk Population, Scenario D: ", my_risky_pop)

        # a list for the pops in risk during this rep:
        time_series.append(my_scenario_D_ts)
        # store: #experiment, max of the series,
        stats.append([i, my_scenario_D_ts.max(),
                      t_range[my_scenario_D_solution.argmax()], my_risky_pop ])
        #del(new_H, my_scenario_D_solution, my_scenario_D_ts,
        #    my_neighbours, my_risky_pop)

    for i in time_series:
        i.plot(color='Crimson', linewidth=3, alpha=.45)
    del(i)

    # I plot the unmodded model for comparison:
    my_model_ts.plot(color='DarkBlue', linewidth=6)
    #pl.title('Scenario D: random Erdos-Renyi networks based on the observed network, red curves.\n Blue curve is the model with the observed connection network ', fontsize=20)
    pl.title('Epidemic scenario three: random Erdos-Renyi networks based on the observed network, red curves.\n Blue curve is the model with the observed connection network ', fontsize=20)
    #pl.title('Scenario D: randomized connections for initial infected island, red curves.'
    #         + r' Blue curve is the model with the observed connection network ', fontsize=20)
    pl.xlabel('time (days)', fontsize=20)
    #pl.legend(loc=1)
    pl.savefig("Erdos-Renyi_scenario_D.png",dpi=400)
    #print(stats)
    #pl.show()

    # histograms of max peak and time to max peak of the experiment:
    histogram_maxatpeak = np.zeros(len(stats), dtype=float)
    histogram_timetomax = np.zeros(len(stats), dtype=float)
    histogram_popatrisk = np.zeros(len(stats), dtype=float)
    for i in range(0, len(stats)):
        histogram_maxatpeak[i] = stats[i][1]
        histogram_timetomax[i] = stats[i][2]
        histogram_popatrisk[i] = np.log(stats[i][3])
        if stats[i][3] == 0:
            histogram_popatrisk[i] = 0 # to correct the problem with -Inf
    del(i)

    pl.figure()
    pl.hist(histogram_maxatpeak, color='crimson', bins=10)
    pl.xlabel('Number of infected at maximum epidemic peak', fontsize=20)
    pl.ylabel('frequency', fontsize=20)
    pl.savefig("infected_at_maximum_epidemic_peak_scenario_D.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist(histogram_timetomax, bins=10)
    pl.xlabel('Time in days to reach epidemic peak', fontsize=20)
    pl.ylabel('frequency', fontsize=20)
    pl.savefig("time_to_peak_scenario_D.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist(histogram_popatrisk, bins=10)
    pl.xlabel('Population at risk on first degree', fontsize=20)
    pl.ylabel('frequency', fontsize=20)
    pl.savefig("riksy_pop_scenario_D.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist2d(histogram_maxatpeak, histogram_timetomax, bins=15)
    pl.xlabel('Number of infected at maximum epidemic peak', fontsize=15)
    pl.ylabel('Time in days to reach epidmic peak', fontsize=15)
    pl.title('2D Histogram', fontsize=22)
    pl.colorbar()
    pl.savefig("Hist1_scenario_D.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist2d(histogram_maxatpeak, histogram_popatrisk, bins=15)
    pl.xlabel('Number of infected at maximum epidemic peak', fontsize=15)
    pl.ylabel('Population at risk on first degree', fontsize=15)
    pl.title('2D Histogram', fontsize=22)
    pl.colorbar()
    pl.savefig("Hist2_scenario_D.png",dpi=400)
    #pl.show()
    pl.figure()
    pl.hist2d(histogram_timetomax, histogram_popatrisk, bins=15)
    pl.xlabel('Time in days to reach epidmic peak', fontsize=15)
    pl.ylabel('Population at risk on first degree', fontsize=15)
    pl.title('2D Histogram', fontsize=22)
    pl.colorbar()
    pl.savefig("Hist3_scenario_D.png",dpi=400)
    #pl.show()

    # 3d-histogram (pending):
    return (sol_matrix)



#################################################################################
#################################################################################
## Code to produce numerical solutions of the network-SEIR infection model     ##
#################################################################################
#################################################################################

## Setting paths:
#setting paths:
if (os.environ['HOME'] == "/home/carlos"):
    data_project_folder='projects/chik-caribbean/'
else:
    data_project_folder='projects/chik-caribbean/'

############### Diseases infection parameters ###################################
# Disease infection parameters:
# mu=1/(70*365.0)  # natality/mortality rate for a demographic model
beta=182.5/365.0 # transmission rate
beta = .78       # (modded)
sigma=1/4.032    # (IIP)
sigma=1/2.032    # (IIP) (modded)
gamma=1/7.194    # recovery rate
gamma=1/4.194    # recovery rate (modded)

beta  = 1.5*beta
sigma = 1.2*sigma
gamma = 1.3*gamma
# further tunning up of the parameters
beta = .69*beta                                                                #
sigma = 1.3*sigma                                                              #
gamma = 1.4*gamma

m_factor = .56e-2
pop_rescale = .001                                                         #

# least squares best fitted parameters for the long observed TS:
beta=0.864832
sigma=0.798620
gamma=0.491582

#least squares best fitted  parameters for the short observed TS:
beta=1.209832
sigma=0.35361999999999999
gamma=0.56810111111111117


#################################################################################
#################################################################################
#################################################################################
# experiment with more realistic parametric initial values condition:
## (done on Sep 17 2017)
beta=0.165
sigma=0.25
gamma=0.14
# optimized parameter values obtained from previous initial conditions:
# (fitted for the long time series)
#beta=0.4561
#sigma=0.385
#gamma=0.2461 
#################################################################################
#################################################################################
#################################################################################



##fitted parameters for the long TS:
#beta=0.448
#sigma=0.655
#gamma=0.269


## least squares best fitted parameters for the long observed TS:
#beta=0.487850
#sigma=0.482278
#gamma=0.280734

# for the least square fitting method:
init_beta = beta
init_sigma = sigma
init_gamma = gamma
#init_beta = round(beta)
#init_sigma = round(sigma)
#init_gamma = round(gamma)
init_bound = 0.7  # for building an interval around the intial value
param_partition = 10.0  # number of points that divide the interval
numberof_recurrences = 1  # number of iterations for the least square fitting

############### network related parameters: #####################################
epi_classes = 5  # SEIRN

############## solver related parameters: #######################################
t_start = 0.0
t_end=45.0*7     # total number of days: number weeks times number of days
#t_end=302   # according obvservation dates 20-10-2013 to 24-08-2014 (302 days)
#t_inc=.10       # time step, time increment
#t_inc=1.0/24.0   # each day between 12 hours
t_inc=1.0   # daily 
t_range = np.arange(t_start, t_end+t_inc, t_inc)


########## experiments related parameters #######################################
reps = 500 # number of experiments or repetition for each scenario
experiment = np.random.random_sample() # naming the exporiment ID
#################################################################################
print ("-------------------------------")
print ("beta = "), beta
print ("sigma = "), sigma
print ("gamma = "), gamma
print ("-------------------------------")
print ("Number of reps: "), reps
print ("-------------------------------")
#################################################################################

#################################################################################
# From 2 datasets generated by Albert's scripts:
# - island_connections_caribbean_3.csv'
# - airports.csv,
# I create a network model of the flights connection in the selected area
# and then I couple a SEIR model with vector-born infection
# The flow is:
# - create a pandas data frame (my_aiports_df) that have ALL aiports in the
# selected area and their connections,
# - create a pandas data frame (my_latlog_df) that have ALL airports in the
# selected area and their coordinates, country to which belongs and a 'political
# region' based roughly in their official spoken language
# - create a networkx graph (my_g) based in the previous objects to produce a network
# model of flight connection of ALL airports,
# - create a pandas dataframe with an adjacency matrix for the he reduced flight
# model BY COUNTRY instead of all airports
#
# Creating a matrix from the df with countries and connections:
#################################################################################

################# FUNCTION CALLS: ###############################################
my_airports_df   = f_airports(data_project_folder) # load my_aiports_df dataframe



my_latlong_df    = f_latlon_df()



my_g             = f_graph(my_airports_df, my_latlong_df)



mydf             = f_dataframe(my_g, my_airports_df, my_latlong_df)



H                = f_reduced_g(mydf, my_g, pop_rescale)



my_observed_chik_ts_indv = f_plot_observed_chik_indv(data_project_folder)



my_observed_chik_ts = f_aggregated_ts_observed_chik(data_project_folder)



# indicate initial infected node for constructing Y0:
node = 'Saint Martin'
#node = 'Colombia' ## testing
Y0 = f_yo_constructor(H, node)


#my_model_plot = f_plot_model(m_factor, H, Y0, t_start, t_end, t_inc, \
#                             t_range, epi_classes, beta, sigma, gamma)
#


#my_model_solution_nodes = f_node_model_solution(m_factor, H, Y0, t_start, t_end,
#                                                t_inc, t_range, epi_classes,
#                                                beta, sigma, gamma)
#
#
#
#my_model_nodes_ts = f_node_model_ts(my_model_solution_nodes)
#my_dict = f_node_model_ts(my_model_solution_nodes)  # do not forget to switch the dict opt in the fnc.
#
#
#
#
#my_model_solution = f_aggregated_model_solution(m_factor, H, Y0, t_start, t_end,
#                                                t_inc, t_range, epi_classes,
#                                                beta, sigma, gamma)
#
#
#
#my_model_ts = f_aggregated_model_ts(my_model_solution)
#
#
#
#my_model_solution_plot = f_model_solution_plot(my_model_solution)



#my_minimum_square_error =  minimum_square_error(my_observed_chik_ts,
#                                                 H, Y0, t_start, t_end, t_inc,
#                                                 t_range, epi_classes,
#                                                 init_beta,init_sigma, init_gamma,
#                                                 init_bound, param_partition, numberof_recurrences, experiment)
#
#
#my_fitted_model = fitted_model(my_observed_chik_ts, H, Y0, t_start, t_end, t_inc, t_range, epi_classes,
#                 init_beta, init_sigma, init_gamma, init_bound,
#                 numberof_recurrences, param_partition)
#

#my_degree_histogram = degree_histogram(H)


## Scenario A (Saint Martin with randomized links):
#my_new_H_scenario_A = f_H_scenario_A(H)
#my_scenario_A_plot = f_plot_scenario_A(reps, H, )
##np.savetxt("scenario_A_sol.csv", my_scenario_A_plot, delimiter=",")
#
#
## Scenario B: randomized seed location, observed network unchanged
#my_randomize_Y0_scenario_B = f_shuffle_Y0_scenario_B(H)
##my_scenario_B_plot = f_plot_scenario_B(reps, H, my_randomize_Y0_scenario_B)
#my_scenario_B_plot = f_plot_scenario_B2(H)
##np.savetxt("scenario_B_sol.csv", my_scenario_B_plot, delimiter=",")


### Scenario C (Saint Martin connected to the rest nodes):
#my_new_H_scenario_C = f_H_scenario_C(H)
#my_scenario_C_plot = f_plot_scenario_C(my_new_H_scenario_C, Y0)
#np.savetxt("scenario_C_sol.csv", my_scenario_C_plot, delimiter=",")


## Scenario D (Erdos-Renyi G(n,p) from the observed flight network):
#my_new_H_scenario_D = f_H_scenario_D(H)
#my_scenario_D_plot = f_plot_scenario_D(reps, my_new_H_scenario_D, Y0)
#np.savetxt("scenario_D_sol.csv", my_scenario_D_plot, delimiter=",")

## writing a report file for the experiment:
#f = open('./figures/experiment=%f' %experiment, 'w')
#f.write('experiment = %f' %experiment)
#f.close()







### WILL COMMENT FROM NOW ON TEMPORARILY ########################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
### ALL THIS NEED TO BE REFACTORED INTO MORE GENERAL FUNCTIONS ##
#################################################################
#################################################################
#################################################################################
#################################################################################
#################################################################################
## Miscelaneous:
## routines to produce the shortened time series of obervations and the solution
## of the model fitted to the shortened observation time series:
#################################################################################
#################################################################################
#################################################################################
#
#
# shorten down the original observed time series:
my_obs_time_series = obs_ts_short(my_observed_chik_ts)

# call the optimizer with the shortened time series as one of its arguments in
# order to fit the model solution to that particular series:

#init_beta =  beta
#init_sigma = sigma
#init_gamma = gamma

init_beta = 0.66
init_sigma = 0.24
init_gamma = 0.14
print("initial values for beta, sigma, gamma", init_beta, init_sigma, init_gamma)
print("")
rescale_factor = 0.7
### Fitting the aggregated model:
#my_fitted_model = fitted_model(my_obs_time_series, H, Y0, t_start, t_end, t_inc, t_range, epi_classes,
#                 init_beta, init_sigma, init_gamma, init_bound,
#                 numberof_recurrences, param_partition)

#################################################################################
# fitting the model using sensible parametric space found in the literature for
# beta, sigma and gamma. beta is unknow as it integrates both mosquito to human
# and human to mosquito infection rates. For this we have a larger uncertainty
# and explore 'blindly' the possible parameter space. For the others parameters
# we build the intervals to explores based in what have be reported by Staples,
# 2009, Pialoux 2007, Boelle 2008 (look at the chik abm paper):
param_partition = 8
rescale_factor = 0.3
rescale_space = np.linspace(0,1,11) # ~(0, 0.1, ..., 1.0)
beta_space = np.linspace(0, 2, param_partition)
sigma_space = np.linspace(0.04667, 0.5, param_partition) # (21.4^-1 days, to, 2^.-1 days))
gamma_space = np.linspace(0.03, 0.2, param_partition) # (33.33^-1 days, to, 5^-1 days))

# Fitting the aggregated model at different scale factors values
# (scale factors multiply the model solution by values from 0,..,1))
#for i in np.arange(0, len(rescale_space)):
#    rescale_factor = rescale_space[i]
#    print()
#    print()
#    print("############################")
#    print("RESCALE FACTOR = ", rescale_factor)
#    print("############################")
#    print()
#    my_fitted_model = fitted_model2(my_obs_time_series, H, Y0, t_start, t_end, t_inc,
#                                    t_range, epi_classes, init_beta, init_sigma, init_gamma,
#                                    init_bound, numberof_recurrences, param_partition)
#
#
#    print("best fitted beta = ", my_fitted_model[0])
#    print("best fitted sigma = ", my_fitted_model[1])
#    print("best fitted gamma = ", my_fitted_model[2])
#    print("")
#    # fitted parameter values:
#    beta = my_fitted_model[0]
#    sigma = my_fitted_model[1]
#    gamma = my_fitted_model[2]
#
#    # Then I call the solver to solve the ODE with the fitted values:
#
#    ## parameter values hand-tuned:
#    #beta = 0.66
#    #sigma = 0.24
#    #gamma = 0.14
#    ##best fitted values according the optimization algorhytm:
#    #beta = 0.46333
#    #sigma = 1.306667
#    #gamma = 0.203333
#    my_model_solution = f_aggregated_model_solution(m_factor, H, Y0, t_start, t_end,
#                                                    t_inc, t_range, epi_classes,
#                                                    beta, sigma, gamma)
#    my_model_ts = f_aggregated_model_ts(my_model_solution)
#
#    # computing the Euclidean distance between my_model_ts and my_obs_time_series.
#    my_model_ts_reduced = f_sample_solution(my_obs_time_series, my_model_ts)
#    distance = f_vector_distance(my_obs_time_series, my_model_ts_reduced)
#
#    # plotting figs:
#    # some configuration first:
#    xlabelsize = 12
#    ylabelsize = 12
#    pl.rcParams['xtick.labelsize'] = xlabelsize
#    pl.rcParams['ytick.labelsize'] = ylabelsize
#    hfont = {'fontname':'Helvetica'}
#    sefont = {'fontname':'Serif'}
#    pl.figure(figsize=(18,7))
#    pl.plot(my_obs_time_series)
#    #pl.plot(my_model_ts*.3)
#    pl.plot(my_model_ts)
#    pl.title('Model fitting -least squares- beta=%f, sigma=%f gamma=%f, dist. between obs and sol=%f, param_partition=%f, rescale factor=%f' %(beta, sigma, gamma, distance, param_partition, rescale_factor))
#    pl.ylabel("Number of infecteds", size=15, **sefont)
#    pl.xlabel("Time",size=15, **sefont)
#    #pl.savefig('./figures/FITTED_AGGREGATED_MODEL_experiment=%f_beta=%f_sigma=%f_gamma=%f_distance=%f_param_partition=%f_rescale_factor=%f.png'
#    #           %(experiment, beta, sigma, gamma, distance, param_partition, rescale_factor), dpi=400)
#    pl.savefig('./figures/FITTED_AGGREGATED_MODEL_experiment=%f_beta=%f_sigma=%f_gamma=%f_distance=%f_param_partition=%f_rescale_factor=%f.pdf'
#               %(experiment, beta, sigma, gamma, distance, param_partition, rescale_factor), dpi=400)
#
#    #pl.figure(figsize=(15,7))
#    #pl.plot(my_obs_time_series)
#    #pl.plot(my_model_ts_reduced)
#    #pl.title('Model fitting (Solution Sampled)-least squares- beta=%f, sigma=%f gamma=%f, dist. between obs and sol=%f' %(beta, sigma, gamma, distance))
#    #pl.ylabel("Number of infecteds", size=15, **sefont)
#    #pl.xlabel("Time",size=15, **sefont)
#    #pl.savefig('./figures/FITTED_AGGREGATED_SAMPLED_MODEL_experiment=%f_beta=%f_sigma=%f_gamma=%f_distance=%f.png'
#    #           %(experiment, beta, sigma, gamma, distance), dpi=400)
#    #pl.savefig('./figures/FITTED_AGGREGATED_SAMPLED_MODEL_experiment=%f_beta=%f_sigma=%f_gamma=%f_distance=%f.pdf'
#    #           %(experiment, beta, sigma, gamma, distance), dpi=400)
#    #pl.show(block=False)
#

################################...###############################################
##################################################################################
## Scripts to get model solution at individual nodes:
## Now I need to fit the parameters to observed incidence in some individual
## nodes:
## vector of initialization:
#node = "Saint Martin"
#Y0_node = np.zeros(5)
#Y0_node[0] = H.node[node]['pop_size'] - 1
#Y0_node[1] = 0
#Y0_node[2] = 1
#Y0_node[3] = 0
#Y0_node[4] = H.node[node]['pop_size']
#
#
## solving the eqs.:
#YSOL_node = spi.odeint(ode_eq, Y0_node, t_range, args=(beta, sigma, gamma))
#YSOL_node = YSOL_node*.15 # only 15% are observed or recorded
#
#YSOL_node_ts = f_aggregated_model_ts(YSOL_node[:,2])
#
## Next I need to fit parameter to an outbreak in a particular node
## Once I have the model solution, YSOL_node_ts,  as a TimeSeries object, I need
## to sample the points of the solution to match those of the observation vector
## so I can measure a (Euclidean) distance between the solution and the
## observation and porceed with the optimization process:
#
## individual observed time series for the node in question:
## fistly, shorten the observed time series of the node of interest such it has
## the same length as the modeled short series:



## Call the functions:
## I take the first 30 data points of the ts as they corresponds to the shorter
## time series analyzed for the first stages of the outbreak
#my_observed_chik_ts = my_observed_chik_ts_indv["Saint-Martin"][0:30]
##my_observed_chik_ts = my_observed_chik_ts_indv["Guadeloupe"][0:30]
##my_observed_chik_ts = my_observed_chik_ts_indv["Martinique"][0:30]
#
#
#my_fitted_model_node = fitted_model_node(my_observed_chik_ts, H, Y0,
#                                         t_start, t_end, t_inc, t_range, epi_classes,
#                                         init_beta, init_sigma, init_gamma, init_bound,
#                                         numberof_recurrences, param_partition, experiment)
#
#
#
################################...###############################################
##################################################################################





#my_model_time_series = model_ts_short(my_model_ts)

# This is the solution sampled at the time of the observations (short time
# series):
#my_model_sampled_time_series = \
#    f_sample_solution(my_obs_time_series, my_model_time_series)
#
#
## obs vs fitted model:
## plotting figs:
## some configuration first:
#xlabelsize = 15
#ylabelsize = 15
#pl.rcParams['xtick.labelsize'] = xlabelsize
#pl.rcParams['ytick.labelsize'] = ylabelsize
#hfont = {'fontname':'Helvetica'}
#sefont = {'fontname':'Serif'}
#
#pl.figure(figsize=(10,5))
#pl.plot(my_obs_time_series, linewidth=4, color="red", alpha=.7)
#pl.plot(my_model_time_series, linewidth=4, color="blue", alpha=.7)
#pl.title("Observed data (red) vs. fitted model", size=20, **sefont)
#pl.ylabel("Chikungunya Cases", size=15, **sefont)
#pl.xlabel("Time in Weeks",size=15, **sefont)
##pl.savefig('./figures/observed_vs_fitted_ST.png', dpi=200)
##pl.savefig('./figures/observed_vs_fitted_ST.jpg', dpi=400)
#
#pl.figure(figsize=(10,5))
#pl.plot(my_obs_time_series, linewidth=4, color="red", alpha=.7)
#pl.title("Observed data", size=20, **sefont)
#pl.ylabel("Chikungunya Cases", size=15, **sefont)
#pl.xlabel("Time in Weeks",size=15, **sefont)
##pl.savefig('./figures/observed_ST.png', dpi=200)
#
#pl.figure(figsize=(10,5))
#pl.plot(my_model_time_series, linewidth=4, color="blue", alpha=.7)
#pl.title("Fitted model", size=20, **sefont)
#pl.ylabel("Chikungunya Cases", size=15, **sefont)
#pl.xlabel("Time in Weeks",size=15, **sefont)
##pl.savefig('./figures/fitted_ST.png', dpi=200)
#
#pl.show(block=False)


#my_plot_obs_and_model = plot_obs_and_model(my_obs_time_series, my_model_time_series)



## Now I need to produce the same time series for both the shortened observation
## and the shortened solution only up to the max epidemic peak, plot them and
## compute the number of infected at peak, the accumulated number of infected at
## the peak and the number of days to the peak:
#my_obs_max = my_obs_time_series.max()
#
## time at peak vec
#time_at_obs_peak  = pd.Series.idxmax(my_obs_time_series)
#
## number of days to peak
#num_days_at_obs_peak = time_at_obs_peak - my_obs_time_series.index.min()
#num_days_at_obs_peak = int(str(num_days_at_obs_peak)[0:3])
#
## accumulated infected at peak (sort of hack):
#accum_infect_at_obs_peak = sum(my_obs_time_series[my_obs_time_series.index.min():pd.Series.idxmax(my_obs_time_series)])
#
##-----------------------------------#
## The same for the model solution:
#my_model_max = my_model_time_series.max()
#
## time at peak vec
#time_at_model_peak  = pd.Series.idxmax(my_model_time_series)
#
## number of days to peak
#num_days_at_model_peak = time_at_model_peak - my_model_time_series.index.min()
#num_days_at_model_peak = int(str(num_days_at_model_peak)[0:3])
#
## accumulated infected at peak (sort of hack):
#accum_infect_at_model_peak = sum(my_model_time_series[my_model_time_series.index.min():pd.Series.idxmax(my_model_time_series)])
#
### Plot code to check and visualize solutions up to epidemic peaks:
#pl.figure()
#my_obs_ts_2peak = (my_obs_time_series[my_obs_time_series.index.min():pd.Series.idxmax(my_obs_time_series)])
#my_obs_ts_2peak.plot(linewidth=3, color="Crimson", alpha=.6)
#
#my_model_ts_2peak = (my_model_time_series[my_model_time_series.index.min():pd.Series.idxmax(my_model_time_series)])
#my_model_ts_2peak.plot(linewidth=3, color="blue", alpha=.6)
#
#pl.title("Visualization of epidemic size up to highest peak, red is observed data, blue fitted model")
#pl.xlabel("Time")
#pl.ylabel("Number of infecteds")
#
#
#################################################################
#################################################################
#################################################################
#################################################################
##### routines for summarizing statistics of the scenarios: #####
#################################################################
#################################################################
#################################################################
##my_scenario_str = "scenario_A_sol.csv"
##my_scenario_str = "scenario_B_sol.csv"
#my_scenario_str = "scenario_C_sol.csv"
##my_scenario_str = "scenario_D_ER_network_sol.csv"
##my_scenario_str = "scenario_D_sol.csv"
##
#scenario =  pd.read_csv(my_scenario_str, sep = ",", header = None)
#scenario = scenario.values # converting df to numpy array
#
## we save the homogeneous scenario as well:
#scenario_homog = pd.read_csv("scenario_C_sol.csv", sep = ",", header = None)
#scenario_homog = scenario_homog.values # converting df to numpy array
#
## convert model scenario to a time series:
#rng = pd.date_range('20-10-2013', periods=scenario.shape[0], freq='H')
#rng_homog = pd.date_range('20-10-2013', periods=scenario_homog.shape[0], freq='H')
#
## create the list:
#scenario_ts = [pd.Series(scenario[:,i], index=rng)
#               for i in np.arange(0, scenario.shape[1])]
#scenario_ts_homog = [pd.Series(scenario_homog[:,i], index=rng_homog)
#               for i in np.arange(0, scenario_homog.shape[1])]
### create list for sampled solutions according observed TS:
##scenario_ts_sampled = [f_sample_solution(my_obs_time_series, scenario_ts[i])
##                       for i in np.arange(0, len(my_obs_time_series))]
## for scenario homogeneous (scenario C) use this code instead:
#scenario_ts_sampled = f_sample_solution(my_obs_time_series, scenario_ts[0])
#
#
## max nuber of infecteds (peak):
#peak_vec =  [max(scenario_ts[i])
#             for i in np.arange(0, len(scenario_ts))]
#peak_vec_homog =  [max(scenario_ts_homog[i])
#                   for i in np.arange(0, len(scenario_ts_homog))]
##peak_vec_sampled_ts =  [max(scenario_ts_sampled[i])
##                            for i in np.arange(0, len(scenario_ts_sampled))]
## for scenario homogeneous (scenario C) use this code instead:
#peak_vec_sampled_ts = max (scenario_ts_sampled)
#
## time at peak vec
#time_at_peak  = [pd.Series.idxmax(scenario_ts[i])
#                 for i in np.arange(0, len(scenario_ts))]
#time_at_peak_homog  = [pd.Series.idxmax(scenario_ts_homog[i])
#                 for i in np.arange(0, len(scenario_ts_homog))]
##time_at_peak_sample_ts  = [pd.Series.idxmax(scenario_ts_sampled[i])
##                 for i in np.arange(0, len(scenario_ts_sampled))]
## for scenario homogeneous (scenario C) use this code instead:
#time_at_peak_sample_ts = pd.Series.idxmax(scenario_ts_sampled)
#
## number of days at peak
#num_days_at_peak = [time_at_peak[i] - scenario_ts[i].index.min()
#                    for i in np.arange(0, len(scenario_ts))]
#num_days_at_peak_homog = [time_at_peak_homog[i] - scenario_ts_homog[i].index.min()
#                    for i in np.arange(0, len(scenario_ts_homog))]
##num_days_at_peak_sampled_ts = [time_at_peak_sample_ts[i] - scenario_ts_sampled[i].index.min()
##                    for i in np.arange(0, len(scenario_ts_sampled))]
## for scenario homogeneous (scenario C) use this code instead:
#num_days_at_peak_sampled_ts = time_at_peak_sample_ts - scenario_ts_sampled.index.min()
#
#
## hack:
#num_days_at_peak = [int(str(num_days_at_peak[i])[0:3])
#                    for i in np.arange(0, len(scenario_ts))]
#num_days_at_peak_homog = [int(str(num_days_at_peak_homog[i])[0:3])
#                    for i in np.arange(0, len(scenario_ts_homog))]
##num_days_at_peak_sampled_ts = [int(str(num_days_at_peak_sampled_ts[i])[0:3])
##                    for i in np.arange(0, len(scenario_ts_sampled))]
## for scenario homogeneous (scenario C) use this code instead:
#num_days_at_peak_sampled_ts = int(str(num_days_at_peak_sampled_ts)[0:3])
#
#
## accumulated infected at peak (sort of hack):
#accum_infect_at_peak = [sum(scenario_ts[i][scenario_ts[i].index.min():pd.Series.idxmax(scenario_ts[i])])
#                        for i in np.arange(0, len(scenario_ts))]
#accum_infect_at_peak_homog = [sum(scenario_ts_homog[i][scenario_ts_homog[i].index.min():pd.Series.idxmax(scenario_ts_homog[i])])
#                        for i in np.arange(0, len(scenario_ts_homog))]
##accum_infect_at_peak_sampled_ts = [sum(scenario_ts_sampled[i][scenario_ts_sampled[i].index.min():pd.Series.idxmax(scenario_ts_sampled[i])])
##                        for i in np.arange(0, len(scenario_ts_sampled))]
## for scenario homogeneous (scenario C) use this code instead:
#accum_infect_at_peak_sampled_ts = sum(scenario_ts_sampled[scenario_ts_sampled.index.min():pd.Series.idxmax(scenario_ts_sampled)])
#
#
#
### Plot code to check and visualize solutions up to eidemic peaks:
#pl.figure(figsize=(12, 7))
#[scenario_ts[i][scenario_ts[i].index.min():pd.Series.idxmax(scenario_ts[i])].plot(linewidth=2, alpha=.15, color='lightcoral')
# for i in np.arange(-1, len(scenario_ts))]
#[scenario_ts_homog[i][scenario_ts_homog[i].index.min():pd.Series.idxmax(scenario_ts_homog[i])].plot(linewidth=6, alpha=.75, color='blue')
# for i in np.arange(-1, len(scenario_ts_homog))]
#
## plot both the observed and the fitted model too:
#my_obs_ts_2peak.plot(linewidth=4, color="Crimson", alpha=.9)
#my_model_ts_2peak.plot(linewidth=4, color="blue", alpha=.7)
#
#pl.title("Visualization of epidemic size up to highest peak, %s"%(my_scenario_str), size=15, **sefont)
#pl.xlabel("Time", size=15, **sefont)
#pl.ylabel("Num. infecteds, pink: scenarios, red: observed and blue: model", size=15, **sefont)
#pl.show(block=False)
#
#
###############################################################################
###############################################################################
###-----------------------------------------#
###-----------------------------------------#
### Choose between plotting histograms 'A' or 'B'
###-----------------------------------------#
###-----------------------------------------#
#### Histograms 'A'
###pl.figure(figsize=(15,9))
#### pl.suptitle("Histograms for %s"%(my_scenario_str))
###pl.subplot(331)
###pl.hist(peak_vec)
###pl.title("Infecteds at highest peak, %s"%(my_scenario_str))
###pl.xlabel("Number of infecteds")
###pl.ylabel("frequency")
###
###pl.subplot(335)
###pl.hist(num_days_at_peak)
###pl.title("Days to highest peak, %s"%(my_scenario_str))
###pl.xlabel("Number of days")
###pl.ylabel("frequency")
###
###pl.subplot(339)
###pl.hist(accum_infect_at_peak)
###pl.title("Inf accumulated up to highest peak, %s"%(my_scenario_str))
###pl.xlabel("Number of (acc.) infecteds")
###pl.ylabel("frequency")
###
###pl.tight_layout()
###pl.savefig("./figures/histograms_A_%s.png"%(my_scenario_str),
###           dpi=400, orientation='landscape')
###
###############################################################################
###############################################################################
#### plotting all together:
#### Histograms 'B'
###pl.figure(figsize=(15,10))
###pl.suptitle("Histograms for %s"%(my_scenario_str))
###pl.subplot(221)
###pl.hist(peak_vec)
###pl.title("Infecteds at highest peak, %s"%(my_scenario_str))
###pl.xlabel("Number of infecteds")
###pl.ylabel("frequency")
###
###pl.subplot(222)
###pl.hist(num_days_at_peak)
###pl.title("Days to highest peak, %s"%(my_scenario_str))
###pl.xlabel("Number of days")
###pl.ylabel("frequency")
###
###pl.subplot(223)
###pl.hist(accum_infect_at_peak)
###pl.title("Inf accumulated up to highest peak, %s"%(my_scenario_str))
###pl.xlabel("Number of (acc.) infecteds")
###pl.ylabel("frequency")
###
###pl.subplot(224)
###[scenario_ts[i][scenario_ts[i].index.min():pd.Series.idxmax(scenario_ts[i])].plot(linewidth=2, alpha=.15, color='lightcoral')
### for i in np.arange(-1, len(scenario_ts))]
#### plot both the observed and the fitted model too:
###my_obs_ts_2peak.plot(linewidth=5, color="Crimson", alpha=.9)
###my_model_ts_2peak.plot(linewidth=5, color="blue", alpha=.7)
###pl.title("Epidemic size up to highest peak (pink: scenarios, red:obs, blue:model), %s"%(my_scenario_str))
###pl.xlabel("Time")
###pl.ylabel("Number of infecteds")
###
###pl.tight_layout()
###pl.savefig("./figures/histograms_B_%s.png"%(my_scenario_str),
###           dpi=400, orientation='landscape')
###pl.savefig("./figures/histograms_B_%s.pdf"%(my_scenario_str),
###           dpi=400, orientation='landscape')
##
##
### for xavi's talk:
##pl.figure(figsize=(12,7))
##[scenario_ts[i][scenario_ts[i].index.min():pd.Series.idxmax(scenario_ts[i])].plot(linewidth=4, alpha=.30, color='lightcoral')
##  for i in np.arange(-1, len(scenario_ts))]
### plot both the observed and the fitted model too:
##my_obs_ts_2peak.plot(linewidth=5, color="Crimson", alpha=.9)
##my_model_ts_2peak.plot(linewidth=5, color="blue", alpha=.7)
##pl.title("Epidemic size up to highest peak (pink: scenarios, red: obs, blue: model), %s"%(my_scenario_str), size=17)
##pl.xlabel("Time", size=15)
##pl.ylabel("Number of infecteds", size=15)
##
###pl.savefig("./figures/scenario%s.png"%(my_scenario_str),
###                      dpi=400, orientation='landscape')
###pl.savefig("./figures/scenario_%s.pdf"%(my_scenario_str),
###           dpi=400, orientation='landscape')
##
##pl.show(block=False)
##
##
###-----------------------------------------#
###-----------------------------------------#
################################################################################
################################################################################
###
### THIS NEEDS TO BE REFACTORED!!!
###
### Miscelaneous code:
## getting the set of nodes belonging to a different language region:
#anglo_nodes = [i for i in H.nodes() if H.node[i]['political_region'] == 'anglo']
#dutch_nodes = [i for i in H.nodes() if H.node[i]['political_region'] == 'dutch']
#latin_nodes = [i for i in H.nodes() if H.node[i]['political_region'] == 'latin']
#french_nodes = [i for i in H.nodes() if H.node[i]['political_region'] == 'french']
## determine what's the node with the larges population of wach set:
#largest_french = max(zip ([ H.node[i]['pop_size'] for i in french_nodes], french_nodes))
#largest_dutch = max(zip ([ H.node[i]['pop_size'] for i in dutch_nodes], dutch_nodes))
#largest_latin = max(zip ([ H.node[i]['pop_size'] for i in latin_nodes], latin_nodes))
#largest_anglo = max(zip ([ H.node[i]['pop_size'] for i in anglo_nodes], anglo_nodes))
#
## Produce a model solution starting in a given location
## in this particular, in a given colored region with the
## largest population
## French:
#node = largest_french[1]
#Y0 = f_yo_constructor(H, node)
#largest_french_sol = f_aggregated_model_solution(m_factor, H, Y0, t_start, t_end,
#                                                t_inc, t_range, epi_classes,
#                                                beta, sigma, gamma)
#
#largest_french_sol_ts = f_aggregated_model_ts(largest_french_sol)
#largest_french_sol_ts_2peak = (largest_french_sol_ts[largest_french_sol_ts.index.min():pd.Series.idxmax(largest_french_sol_ts)])
#
## Ducth
#node = largest_dutch[1]
#Y0 = f_yo_constructor(H, node)
#largest_dutch_sol = f_aggregated_model_solution(m_factor, H, Y0, t_start, t_end,
#                                                t_inc, t_range, epi_classes,
#                                                beta, sigma, gamma)
#
#largest_dutch_sol_ts = f_aggregated_model_ts(largest_dutch_sol)
#largest_dutch_sol_ts_2peak = (largest_dutch_sol_ts[largest_dutch_sol_ts.index.min():pd.Series.idxmax(largest_dutch_sol_ts)])
#
## Latin
#node = largest_latin[1]
#Y0 = f_yo_constructor(H, node)
#largest_latin_sol = f_aggregated_model_solution(m_factor, H, Y0, t_start, t_end,
#                                                t_inc, t_range, epi_classes,
#                                                beta, sigma, gamma)
#
#largest_latin_sol_ts = f_aggregated_model_ts(largest_latin_sol)
#largest_latin_sol_ts_2peak = (largest_latin_sol_ts[largest_latin_sol_ts.index.min():pd.Series.idxmax(largest_latin_sol_ts)])
#
#
## Anglo
#node = largest_anglo[1]
#Y0 = f_yo_constructor(H, node)
#largest_anglo_sol = f_aggregated_model_solution(m_factor, H, Y0, t_start, t_end,
#                                                t_inc, t_range, epi_classes,
#                                                beta, sigma, gamma)
#
#largest_anglo_sol_ts = f_aggregated_model_ts(largest_anglo_sol)
#largest_anglo_sol_ts_2peak = (largest_anglo_sol_ts[largest_anglo_sol_ts.index.min():pd.Series.idxmax(largest_anglo_sol_ts)])
#
## Saint Martin
#node = "Saint Martin"
#Y0 = f_yo_constructor(H, node)
#saint_martin_sol = f_aggregated_model_solution(m_factor, H, Y0, t_start, t_end,
#                                                t_inc, t_range, epi_classes,
#                                                beta, sigma, gamma)
#
#saint_martin_sol_ts = f_aggregated_model_ts(saint_martin_sol)
#saint_martin_sol_ts_2peak = (saint_martin_sol_ts[saint_martin_sol_ts.index.min():pd.Series.idxmax(saint_martin_sol_ts)])
#
#
## plot
#pl.figure(figsize=(10, 5))
#pl.plot(largest_french_sol_ts, linewidth=3, color="blue", alpha=.6)
#pl.plot(largest_dutch_sol_ts, linewidth=3, color="orange", alpha=.6)
#pl.plot(largest_latin_sol_ts, linewidth=3, color="brown", alpha=.6)
#pl.plot(largest_anglo_sol_ts, linewidth=3, color="green", alpha=.6)
#pl.plot(saint_martin_sol_ts, linewidth=6, color="black")
#pl.title("Black=Saint Martin, Blue=Guadaloupe, Orange=Curacao, Green=Jamaica, Brown=Colombia")
#pl.ylabel("Chik Cases")
#pl.show(block=False)
#
## saving the ts
#largest_french_sol_ts.to_csv('largest_french_ts.csv')
#largest_dutch_sol_ts.to_csv('largest_dutch_ts.csv')
#largest_latin_sol_ts.to_csv('largest_latin_ts.csv')
#largest_anglo_sol_ts.to_csv('largest_anglo_ts.csv')
#
#
### plot paper's figure 3
#[scenario_ts[i][scenario_ts[i].index.min():pd.Series.idxmax(scenario_ts[i])].plot(linewidth=3, alpha=.35, color='lightcoral' )
# for i in np.arange(-1, len(scenario_ts))]
#[scenario_ts_homog[i][scenario_ts_homog[i].index.min():pd.Series.idxmax(scenario_ts_homog[i])].plot(linewidth=6, alpha=.75, color='blue')
# for i in np.arange(-1, len(scenario_ts_homog))]
#saint_martin_sol_ts_2peak.plot(linewidth=6, color="black")
## little trick to have both color and linetype
#largest_french_sol_ts_2peak.plot(linewidth=5, color="blue")
#largest_dutch_sol_ts_2peak.plot(linewidth=5, color="orange")
#largest_latin_sol_ts_2peak.plot(linewidth=5, color="brown")
#largest_anglo_sol_ts_2peak.plot(linewidth=5, color="green")
##
#largest_french_sol_ts_2peak.plot(linewidth=2, linestyle='-.', color="black")
#largest_dutch_sol_ts_2peak.plot(linewidth=2, linestyle='-.', color="black")
#largest_latin_sol_ts_2peak.plot(linewidth=2, linestyle='-.', color="black")
#largest_anglo_sol_ts_2peak.plot(linewidth=2, linestyle='-.', color="black")
#saint_martin_sol_ts_2peak.plot(linewidth=6, color="black")
## plotting figs:
## some configuration first:
#xlabelsize = 15
#ylabelsize = 15
#pl.rcParams['xtick.labelsize'] = xlabelsize
#pl.rcParams['ytick.labelsize'] = ylabelsize
#hfont = {'fontname':'Helvetica'}
#sefont = {'fontname':'Serif'}
#pl.ylabel("Number of infecteds", size=15, **sefont)
##pl.xlabel("Time",size=15, **sefont)
#pl.show(block=False)
#
#
################################################################################
########################## INTERVENTIONS #######################################
### Code to plot intervention solutions up to their global peaks
## 1. load the solution list of the intervenation scenarios:
## the list of solutions of interventions in stored in
## ./intervention_scenarios_05_16_2017/my_solution_list.txt
## to read the saved list:
#import pickle
#with open("./intervention_scenarios_05_16_2017/my_solution_list.txt", "rb") as fp:
#    solution_list = pickle.load(fp)
#with open("./intervention_scenarios_05_16_2017/my_time_vec.txt", "rb") as fp:
#    my_time_vec = pickle.load(fp)
#
## plot
#MAP = 'winter'
##MAP = 'RdYlBu'
#cmap = pl.get_cmap(MAP)
#color = [cmap(i) for i in np.linspace(0, 1, 11)]
#interv_days_list = list([5, 10, 20, 40, 60, 80, 100, 120, 140, 150])
#
## plotting figs:
## some configuration first:
#xlabelsize = 15
#ylabelsize = 15
#pl.rcParams['xtick.labelsize'] = xlabelsize
#pl.rcParams['ytick.labelsize'] = ylabelsize
#hfont = {'fontname':'Helvetica'}
#sefont = {'fontname':'Serif'}
#interv_days_list = list([5, 10, 20, 40, 60, 80, 100, 120, 140, 150])
#
#for i in range(0, len(solution_list)):
#    pl.figure(i, figsize=(10,5))
#    for j in range(0, solution_list[i].shape[1]):
#        max_idx = solution_list[i][:,j].argmax()
#        pl.plot(my_time_vec[0:max_idx],
#                solution_list[i][0:max_idx,j], color=color[j], lw=3, alpha=.7)
#        pl.ylabel("Chikungunya Cases", size=15, **sefont)
#        pl.xlabel("Days", size=15, **sefont)
#        pl.tight_layout()
#    pl.axvline(x = interv_days_list[i], color="r", linestyle="--", lw=3)
#    pl.show()
#
#
### This piece of code dosn't work. This is workk pending/in progress
##solution_list_ts = []
##solution_list_ts_2peak = []
##for i in range(0, len(solution_list)): # loop over the intervention days (10)
##    ts = list()
##    ts_2peak = list()
##    # loop over the intervention levels (11 levels)
##    for j in range(0, solution_list[i].shape[1]):
##        # print(f_aggregated_model_ts(solution_list[i][:,j]))
##        solution_list[i][:,j] = f_aggregated_model_ts(solution_list[i][:,j])
##        ts.append(pd.core.series.Series(f_aggregated_model_ts(solution_list[i][:,j])))
##        ts_2peak.append(ts[j][ts[j].index.min():pd.Series.idxmax(ts[j])])
##
##    solution_list_ts.append(ts)
##    solution_list_ts_2peak.append(ts_2peak)
##
### plot
##MAP = 'winter'
###MAP = 'RdYlBu'
##cmap = pl.get_cmap(MAP)
##color = [cmap(i) for i in np.linspace(0, 1, len(solution_list_ts[0])+1)]
##
### plotting figs:
### some configuration first:
##xlabelsize = 15
##ylabelsize = 15
##pl.rcParams['xtick.labelsize'] = xlabelsize
##pl.rcParams['ytick.labelsize'] = ylabelsize
##hfont = {'fontname':'Helvetica'}
##sefont = {'fontname':'Serif'}
### changing the list interv_days_list to pd.series object
### interv_days_list = list([5, 10, 20, 40, 60, 80, 100, 120, 140, 150])
##interv_days_list_date = list(["2013-10-25",
##                              "2013-11-09",
##                              "2013-10-30",
##                              "2013-11-29",
##                              "2013-12-19",
##                              "2014-01-08",
##                              "2014-01-28",
##                              "2014-02-17",
##                              "2014-03-09",
##                              "2014-03-03"])
##interv_days_list_series = list()
##for i in  range (0, len(interv_days_list_date)):
##    interv_days_list_series.append(pd.to_datetime(interv_days_list_date[i]))
##
##
### loop over the intervention days (10 days)
##for i in range(0, len(solution_list_ts_2peak)):
##    pl.figure(i, figsize=(10,5))
##    # loop over intervention levels of day j (11 levels):
##    for j in range(0, len(solution_list_ts_2peak[i])):
##        solution_list_ts_2peak[i][j].plot(color=color[j], linewidth=4, alpha=.7)
##        pl.axvline(x = interv_days_list_series[i], color="r", linestyle="--", lw=2)
##        pl.ylabel("Chikungunya Cases", size=15, **sefont)
##        pl.tight_layout()
##        pl.show(block=False)
##pl.show()

################################################################################
################################################################################




################################################################################
### code for some special plots:
### The uni-dimensional wave:

## build a list with node names and time of the largets peak:

#my_model_nodes_ts # list of 30 elements with ts of infecteds
#
## try another idea with dicts:
#my_dict = f_node_model_ts(my_model_solution_nodes)  # do not forget to switch the dict opt in the fnc.
#my_dict_node_time_and_peak = {}
#for i in my_dict: # my_dict is a dictionary with nodename and infected ts
#    d = {}
#    nodename = i
#    my_dict_node_time_and_peak[i] = d
#    d["peak"] =  my_dict[i].max()
#    d["time"] = 2
#    d["norm_peak"] = 3
#    #if i == "Saint Martin":
#
#
## Code to construct an object that contains both nodes names and thee time where
## maximum occurs for both observation and modeled time series
#
#my_node = []        # nodes name
#my_max = []         # max of the ts
#my_maxtime = []     # time when max occurs
#
#node_model = [i for i in my_model_nodes_ts]
#node_obs = [i for i in my_observed_chik_ts_indv]
#max_model = [my_model_nodes_ts[i].max() for i in my_model_nodes_ts]
#max_obs = [my_observed_chik_ts_indv[i].max() for i in my_observed_chik_ts_indv]
#maxtime_model = [my_model_nodes_ts[i].idxmax() for i in my_model_nodes_ts]
#maxtime_obs = [my_observed_chik_ts_indv[i].idxmax() for i in my_observed_chik_ts_indv]
#
#pairwise_model_max = zip(node_model, maxtime_model)
#
#pairwise_obs_max = zip(node_obs, maxtime_obs)
#
## now I convert the lists to data frames:
#pairwise_model_max_df = pd.DataFrame(pairwise_model_max, columns=["node", "time"])
#pairwise_obs_max_df = pd.DataFrame(pairwise_obs_max, columns=["node", "time"])
#
## sort out the DFs by time
#pairwise_model_max_df.sort_values(by = 'time')
#pairwise_obs_max_df.sort_values(by = 'time')
#
## plotting the time lines for local outbreak maxima:
#pl.figure(figsize=(18,9))
#pl.plot(my_model_time_series, color="red", linewidth=4)
#for i in range(1, pairwise_model_max_df.shape[0]):
#    print( '{0}. {1}'.format(i, str(pairwise_model_max_df.iloc[i]["node"])))
#    # vertical line
#    pl.axvline(x = str(pairwise_model_max_df.iloc[i]["time"]), color = "blue",
#               label = '{0}. {1}'.format(i, str(pairwise_model_max_df.iloc[i]["node"])), alpha=.5)
#    #pl.annotate( str(pairwise_model_max_df.iloc[i]["node"]),
#    #            (str(pairwise_model_max_df.iloc[i]["time"]), np.random.randint(0,400)),
#    #            xytext=(15, 15), textcoords='offset points', arrowprops=dict(arrowstyle='-|>') )
#
#    pl.annotate( i,
#                (str(pairwise_model_max_df.iloc[i]["time"]), np.random.randint(0,450)),
#                xytext=(15, 15), textcoords='offset points', arrowprops=dict(arrowstyle='-|>') )
#
#pl.ylabel("Cases", fontsize=15)
#pl.legend(loc = 2)
#pl.tight_layout()
#pl.show(block = False)
## please save this figure manually rather than programatically as a pdf
## and the to get the best layout for the legends to be seens press Ctrl + F
## before saving it.
#
end = time.time()
print ("----------")
print ("execution time: ", end - start, " secs")
print ("converted to minutes: ", (end - start)/60., "minutes")
print ("converted to hours: ", ((end - start)/60.)/60., "hours")
print ("----------")
del start, end
