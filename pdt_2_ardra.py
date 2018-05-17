# Python script to load data from PDT cross section files 

# Import python packages
import numpy as np
import sys, os, time, inspect

###############################################################################
def main(pdt_cx_list, zaid_list, silo_name):
    # Reads command line with all PDT cx files that need to be converted to SILO. 
    # Example: "python load_pdt_xs.py 'u235.cx u238.cx'  '92235 92238'" 
    # would create the list pdt_cx_list=['u235.cx', 'u238.cx']  
    # and the list zaid_list=['92235', '92238']
    print
    print pdt_cx_list

    # Read PDT cross section files
    isotope = []
    for i in range(len(pdt_cx_list)): 
        cx = cxdata()
        cx.cxfilename = pdt_cx_list[i]
        cx.zaid = zaid_list[i]
        cx.getdata(pdt_cx_list[i])
        isotope.append(cx)

    num_isotopes = len(isotope)

    for i in range(num_isotopes):
        print isotope[i].cxfilename   
        print("    cx.num_ngrps = %i" %isotope[i].num_ngrps )   
        print("    cx.nlgndr = %i" %isotope[i].nlgndr)

        #print(isotope[i].emid)

    # Directory name containing all cross section data
    dir_name = 'pdt_cx_dir'
    os.system('rm -rf %s' %dir_name)
    os.system('mkdir %s' %dir_name)

    fid =open(os.path.join(dir_name, 'num_ngrps'), 'w')
    fid.write('1\n')
    fid.write('%i' %isotope[0].num_ngrps)
    fid.close()

    fid =open(os.path.join(dir_name, 'num_ggrps'), 'w')
    fid.write('1\n')
    fid.write('0')
    fid.close()

    fid =open(os.path.join(dir_name, 'num_grps'), 'w')
    fid.write('1\n')
    fid.write('%i' %isotope[0].num_ngrps)
    fid.close()

    fid =open(os.path.join(dir_name, 'nlgndr'), 'w')
    fid.write('%i\n' %isotope[0].nlgndr)
    for leg_moment in range(isotope[0].nlgndr):
        fid.write('  1\n')    
    fid.close()

    fid =open(os.path.join(dir_name, 'emid'), 'w')
    fid.write('%s\n' %isotope[0].num_ngrps)
    for g in range(isotope[0].num_ngrps):
        fid.write('%20.12e\n' %isotope[0].emid[g])
    fid.close()   

    fid =open(os.path.join(dir_name, 'spgrp'), 'w')
    fid.write('%s\n' %isotope[0].num_ngrps)
    for g in range(isotope[0].num_ngrps):
        fid.write('%20.12e\n' %isotope[0].spgrp[g])
    fid.close()   

    fid =open(os.path.join(dir_name, 'de'), 'w')
    fid.write('%s\n' %isotope[0].num_ngrps)
    for g in range(isotope[0].num_ngrps):
        fid.write('%20.12e\n' %isotope[0].de[g])
    fid.close() 

    fid =open( os.path.join('iso_file'), 'w')
    fid.write('%i\n' %num_isotopes)
    for i in range(num_isotopes): 
        fid.write('/%s\n' %isotope[i].zaid)
    fid.close() 

    for i in range(num_isotopes):  
        isotope_dir = os.path.join(dir_name, isotope[i].zaid)
        os.system('mkdir '+isotope_dir)

        #fid =open( os.path.join(isotope_dir, 'atwt'), 'w')
        #fid.write('1\n')
        #atwt = get_mass(isotope[i].zaid)
        #fid.write('%20.12e\n' %atwt )
        #fid.close() 

        os.system('ln -s atwt '+os.path.join(isotope_dir, 'mass'))

        fid =open( os.path.join(isotope_dir, 'sigt'), 'w')
        fid.write('%i\n' %isotope[i].num_ngrps)
        for g in range(isotope[i].num_ngrps):
            fid.write('%20.12e\n' %isotope[i].sigt[g])
        fid.close() 

        # Ask Peter Brown why siga
        fid =open( os.path.join(isotope_dir, 'siga'), 'w')
        fid.write('%i\n' %isotope[i].num_ngrps)
        for g in range(isotope[i].num_ngrps):
            fid.write('%20.12e\n' %isotope[i].siga[g])
        fid.close() 

        #fid =open( os.path.join(isotope_dir, 'chi'), 'w')
        #fid.write('%i\n' %isotope[i].num_ngrps)
        #for g in range(isotope[i].num_ngrps):
        #    fid.write('%20.12e\n' %isotope[i].chi[g])
        #fid.close() 

        #fid =open( os.path.join(isotope_dir, 'sigfnu'), 'w')
        #fid.write('%i\n' %isotope[i].num_ngrps)
        #for g in range(isotope[i].num_ngrps):
        #    fid.write('%20.12e\n' %isotope[i].sigfnu[g])
        #fid.close() 

        fid =open( os.path.join(isotope_dir, 'n2n'), 'w')
        fid.write('%i\n' %isotope[i].num_ngrps)
        for g in range(isotope[i].num_ngrps):
            fid.write('0.\n')
        fid.close() 

        fid =open( os.path.join(isotope_dir, 'n3n'), 'w')
        fid.write('%s\n' %isotope[i].num_ngrps)
        for g in range(isotope[i].num_ngrps):
            fid.write('0.\n')
        fid.close() 

        fid =open( os.path.join(isotope_dir, 'n4n'), 'w')
        fid.write('%i\n' %isotope[i].num_ngrps)
        for g in range(isotope[i].num_ngrps):
            fid.write('0.\n')
        fid.close()

        for leg_moment in range(isotope[i].nlgndr):
            print leg_moment

            fid =open( os.path.join(isotope_dir, 'p'+str(leg_moment)), 'w')
            fid.write('%i\n' %isotope[i].num_ngrps **2)
            for gp in range(isotope[i].num_ngrps):
                for g in range(isotope[i].num_ngrps):
                    fid.write('%20.12e\n' %isotope[i].sigs_gp_to_g[gp,g,leg_moment])
            fid.close() 

        fid =open( os.path.join(isotope_dir, 'sigf'), 'w')
        fid.write('%i\n' %isotope[i].num_ngrps **2)
        for gp in range(isotope[i].num_ngrps):
            for g in range(isotope[i].num_ngrps):
                fid.write('%20.12e\n' %isotope[i].sigf[gp,g])
        fid.close() 

        if len(isotope[i].decay_const):

            fid =open( os.path.join(isotope_dir, 'num_precursors'), 'w')
            fid.write('1\n')
            fid.write('  %i\n' %isotope[i].num_precursors)
            fid.close()  

            fid =open( os.path.join(isotope_dir, 'lambda'), 'w')
            fid.write('%i\n' %isotope[i].num_precursors)
            for f in range(isotope[i].num_precursors):
                fid.write('%20.12e\n' %isotope[i].decay_const[f])
            fid.close() 

            fid = open( os.path.join(isotope_dir, 'chid'), 'w')
            print isotope[i].num_precursors, isotope[i].num_ngrps
            fid.write('%i\n' %int(isotope[i].num_precursors*isotope[i].num_ngrps) )
            for f in range(isotope[i].num_precursors):
                for g in range(isotope[i].num_ngrps):
                    fid.write('%20.12e\n' %isotope[i].chid[f,g])
            fid.close() 
            print "sum chid =", np.sum(isotope[i].chid)

            fid = open( os.path.join(isotope_dir, 'beta_fnubar'), 'w')
            print isotope[i].num_precursors, isotope[i].num_ngrps
            fid.write('%i\n' %int(isotope[i].num_precursors*isotope[i].num_ngrps) )
            for f in range(isotope[i].num_precursors):
                for g in range(isotope[i].num_ngrps):
                    fid.write('%20.12e\n' %isotope[i].beta_fnubar[f,g])
            fid.close() 
            print "sum beta_fnubar =",  np.sum(isotope[i].beta_fnubar)



    prompt = 0 # Use total nubar
    #prompt = 1 # Use prompt nubar, and include delayed neutron data

    plot_type=0 # Don't draw any plots
    ilog=0 # Only used for plots


    # Silo file name containing all cross section data from the dir_name
    # directory 
    os.system('rm -rf %s' %silo_name)

    # Write input file readme for make_silo_file program
    fid =open('readme_file','w')
    fid.write('Cross sections generated using NJOY with FEDS\n')
    fid.close()

    # Run make_silo_file
    #currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    #os.system(currentdir+'/make_silo_file %s %s readme_file iso_file' %(dir_name,silo_name) )

###############################################################################





###############################################################################
def readfile(filename):
    # Reads a file as a 1D and 2D list 
    # The 1D list contains lines the file, the 2D list delimits each line by space
    fid =open(filename)
    list_1D = []
    list_2D = []
    for line in fid:
        list_1D.append(line)
        list_2D.append(line.split())
    fid.close()
    return list_1D, list_2D
###############################################################################





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class cxdata():
    # Stores data from PDT files as a python class
    def __init__(self):
        self = self

    def getdata(self, cxfile):
        (self.cxfile_list, self.cxfile_array) = readfile(cxfile)

        self.grp_bounds = []
        self.sigt = []
        #self.sigs_elastic = []
        #self.sigs_inelastic = []
        self.sigf_tot = []
        self.inv_spgrp = []
        self.chi = []
        self.chid = []
        self.sigfnu = []
        self.decay_const = []
        self.beta_fnubar = []

        if os.path.isfile(cxfile) == True:
            # Get group structure from within cxfile
            for i in range(len(self.cxfile_list)):           
                if "Scattering order" in self.cxfile_list[i]:
                    self.nlgndr = int(self.cxfile_array[i][2]) + 1

                elif "Group boundaries in eV:" in self.cxfile_list[i]:
                    while 'T =' not in self.cxfile_list[i+1]:
                        i += 1 
                        for j in range(len(self.cxfile_array[i])):
                            self.grp_bounds.append( np.float(self.cxfile_array[i][j])*1e-6 )

                    self.num_ngrps = len(self.grp_bounds) - 1
                    self.num_cxfile_rows = int( np.ceil( np.float(self.num_ngrps) / 5. ) )
                    self.sigs = np.zeros(( self.num_ngrps )) 
                    self.sigs_gp_to_g = np.zeros(( self.num_ngrps, self.num_ngrps, self.nlgndr )) 
                    #self.sigf_tot = np.zeros(( self.num_ngrps ))
                    self.sigf = np.zeros(( self.num_ngrps, self.num_ngrps ))

                elif "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "1":
                    last_line = i + self.num_cxfile_rows
                    #print self.num_cxfile_rows, last_line
                    while i < last_line:
                        i += 1
                        #print i
                        #print self.cxfile_list[i]
                        for j in range(len(self.cxfile_array[i])):                       
                            self.sigt.append( np.float(self.cxfile_array[i][j]) )

                #elif "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "2":
                #    last_line = i + self.num_cxfile_rows 
                #    while i < last_line:
                #        i += 1
                #        for j in range(len(self.cxfile_array[i])):                       
                #            self.sigs_elastic.append( np.float(self.cxfile_array[i][j]) )

                #elif "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "3":
                #    last_line = i + self.num_cxfile_rows
                #    while i < last_line:
                #        i += 1
                #        for j in range(len(self.cxfile_array[i])):                       
                #            self.sigs_inelastic.append( np.float(self.cxfile_array[i][j]) ) 

                elif "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "18":
                    last_line = i + self.num_cxfile_rows
                    while i < last_line:
                        i += 1
                        for j in range(len(self.cxfile_array[i])):                       
                            self.sigf_tot.append( np.float(self.cxfile_array[i][j]) )

                elif "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "259":
                    last_line = i + self.num_cxfile_rows
                    while i < last_line:
                        i += 1
                        for j in range(len(self.cxfile_array[i])):                       
                            self.inv_spgrp.append( np.float(self.cxfile_array[i][j]) )

                elif "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "1018":
                    last_line = i + self.num_cxfile_rows
                    while i < last_line:
                        i += 1
                        for j in range(len(self.cxfile_array[i])):                       
                            self.chi.append( np.float(self.cxfile_array[i][j]) )

                elif "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "1452":
                    last_line = i + self.num_cxfile_rows
                    while i < last_line:
                        i += 1
                        for j in range(len(self.cxfile_array[i])):                       
                            self.sigfnu.append( np.float(self.cxfile_array[i][j]) )      

                #elif "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "2501,":
                elif "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "2519,":
                    leg_moment = int(self.cxfile_array[i][3])
                    while len(self.cxfile_array[i+1]) > 0. and "MT" not in self.cxfile_list[i+1]:
                        i += 1
                        if self.cxfile_array[i][0] == "Sink,":
                            g = int( self.cxfile_array[i][3] )
                            gp = int( self.cxfile_array[i][4] )

                            while len(self.cxfile_array[i+1]) > 0. and self.cxfile_array[i+1][0] not in ["Sink,", "MT"]:
                                i += 1
                                for j in range(len(self.cxfile_array[i])):                  
                                    self.sigs_gp_to_g[gp, g, leg_moment] = np.float(self.cxfile_array[i][j])      
                                    if leg_moment == 0: self.sigs[gp] += np.float(self.cxfile_array[i][j]) 
                                    gp += 1   

                elif "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "1054":
                    i += 1
                    self.num_precursors = int(self.cxfile_array[i][5])

                    while len(self.cxfile_array[i+1]) > 0. and self.cxfile_array[i+1][0] != "MT":
                        i += 1
                        for j in range(len(self.cxfile_array[i])):               
                            self.decay_const.append( np.float(self.cxfile_array[i][j]) * 1e-6 )     

                elif "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "455":
                    last_line = i + self.num_cxfile_rows
                    while i < last_line:
                        i += 1
                        for j in range(len(self.cxfile_array[i])):                       
                            self.beta_fnubar.append( np.float(self.cxfile_array[i][j]) )
 
                elif "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "2055":
                    i += 1
                    f = 0
                    print cxfile     
                    print " i=", i, "f=", f, "array=", self.cxfile_array[i]
                    self.num_precursors = int(self.cxfile_array[i][5])
                    self.chid = np.zeros(( self.num_precursors, self.num_ngrps ))
                    while len(self.cxfile_array[i+1]) > 0. and "MT" not in self.cxfile_list[i+1]:
                        i += 1
                        print "> i=", i, "f=", f, "array=", self.cxfile_array[i]
                        if self.cxfile_array[i][0] == "DNG":
                            f = int(self.cxfile_array[i][1])
                            g = 0

                            while len(self.cxfile_array[i+1]) > 0. and self.cxfile_array[i+1][0] not in ["DNG", "MT"]:
                                i += 1
                                print ">> i=", i, "f=", f, ", g=", g, "array=", self.cxfile_array[i]
                                for j in range(len(self.cxfile_array[i])):                  
                                    self.chid[f, g] = np.float(self.cxfile_array[i][j])      
                                    g += 1  
                

                elif "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "2518,":
                    while len(self.cxfile_array[i+1]) > 0. and "MT" not in self.cxfile_list[i+1]:
                        i += 1
                        if self.cxfile_array[i][0] == "Sink,":
                            g = int( self.cxfile_array[i][3] )
                            gp = int( self.cxfile_array[i][4] )

                            while len(self.cxfile_array[i+1]) > 0. and self.cxfile_array[i+1][0] not in ["Sink,", "MT"]:
                                i += 1
                                for j in range(len(self.cxfile_array[i])):                  
                                    self.sigf[gp, g] = np.float(self.cxfile_array[i][j])      
                                    gp += 1   




        if os.path.isfile(cxfile) == True and len(self.sigfnu) <= 1.:
            # if MT 1452 is not in the PDT file, MT 18 and MT 452 will be combined
            for i in range(len(self.cxfile_list)):  
                if "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "18":
                    last_line = i + self.num_cxfile_rows
                    while i < last_line:
                        i += 1
                        for j in range(len(self.cxfile_array[i])):   
                            # adding 'sigf' value to the 'sigfnu' list
                            # note that, we'll multiply by 'nu' later                  
                            self.sigfnu.append( np.float(self.cxfile_array[i][j]) )

                elif "MT" in self.cxfile_list[i] and self.cxfile_array[i][1] == "452":
                    last_line = i + self.num_cxfile_rows
                    while i < last_line:
                        i += 1
                        for j in range(len(self.cxfile_array[i])):      
                            # multiplying 'sigf' by 'nu' to get 'sigfnu'                 
                            self.sigfnu *= np.float(self.cxfile_array[i][j])

        self.grp_bounds = np.array(self.grp_bounds)
        self.sigt = np.array(self.sigt)
        self.inv_spgrp = np.array(self.inv_spgrp)
        #self.sigs_elastic = np.array(self.sigs_elastic)
        #self.sigs_inelastic = np.array(self.sigs_inelastic)
        self.chi = np.array(self.chi)
        self.sigfnu = np.array(self.sigfnu)
        self.decay_const = np.array(self.decay_const)
        self.sigf_tot = np.array(self.sigf_tot) 
        self.beta_fnubar = np.prod([np.array(self.beta_fnubar), self.sigf_tot], axis=0)

        if len(self.decay_const) != 0:
            new_beta_fnubar = np.zeros((self.num_precursors, self.num_ngrps))
            for f in range(self.num_precursors):
                sum_chid_f = np.sum(self.chid[f])
                self.chid[f] /= sum_chid_f
                new_beta_fnubar[f] = self.beta_fnubar * sum_chid_f
            self.beta_fnubar = new_beta_fnubar          

        #if len(self.sigs_elastic) == 0:
        #    self.sigs_elastic = np.zeros( (self.num_ngrps) )
        #elif len(self.sigs_inelastic) == 0:
        #    self.sigs_inelastic = np.zeros( (self.num_ngrps) )
        #self.sigs = self.sigs_elastic + self.sigs_inelastic

        self.siga = self.sigt - self.sigs

        self.emid  = np.zeros(self.num_ngrps)
        self.de    = np.zeros(self.num_ngrps)    
        #self.spgrp = np.zeros(self.num_ngrps)
        for i in range(self.num_ngrps):
            self.emid[i] = 0.5*(self.grp_bounds[i] + self.grp_bounds[i+1])
            self.de[i]   = self.grp_bounds[i] - self.grp_bounds[i+1]
            #self.spgrp[i] = np.sqrt( 2.*self.emid[i] / 939.57 ) * 2.9979e10 #speed in [cm/s]
            #self.spgrp[i] = np.sqrt( 2.*self.emid[i] / 939.57 ) * 2.9979e4 #speed in [cm/microsecond] (ask Peter Brown)
            
        self.spgrp = 1e-4/self.inv_spgrp

 
        return self
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


###############################################################################
def get_mass(zaid):
    (bdfls_list, bdfls_array) = readfile('/usr/gapps/data/nuclear/endl_official/endl2009.2/bdfls')

    mass = 0.
    for j in range(len(bdfls_list)):
        if "nuclear masses" in bdfls_list[j]:
            for k in range(j, len(bdfls_list)):
                if bdfls_array[k][0] == str(zaid):
                    mass = np.float( bdfls_array[k][1] )
                    print zaid, mass
                    break

    return mass
###############################################################################


#pdt_cx_list = sys.argv[1].split() 
#zaid_list = sys.argv[2].split()

pdt_cx_list = ['Natural_Iron_11t_25r_100b_100u.cx', 'Pure_U235_11t_25r_100b_100u.cx']
zaid_list = ['26000', '92235']

silo_name = 'xs_FEDS_iron_400g_1600e.silo'

main(pdt_cx_list, zaid_list, silo_name)

