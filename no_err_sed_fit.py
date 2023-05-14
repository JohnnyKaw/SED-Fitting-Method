import numpy as np
import astroquery
from astroquery.sdss import SDSS
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs
from astroquery.ipac.irsa import Irsa
from astroquery.ipac.irsa.irsa_dust import IrsaDust
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy import coordinates as coords
import astropy.units as u
import matplotlib.pyplot as plt
from numpy import random
import csv
import scipy
import os
import sys
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from numpy import loadtxt
from numpy import random
from numpy import sqrt
from numpy import log
from numpy import array
from numpy import sum
from numpy import exp
from scipy.interpolate import interp1d
import pyphot
from pyphot.phot import unit
from pyphot.phot import Filter
from multiprocessing import Pool


			#read in filter transmissivity and create filters for synthetic photometry
os.chdir('/media/intern2/DD82-B92F/RESEARCH_FALL_2022/FILTERS/trimmed_transmissivity_curves/')
#gaia_g_lams = np.loadtxt('gaiaG.TRM')[:,0] * unit['AA']
#gaia_g_trans = np.loadtxt('gaiaG.TRM')[:,1]
J_2MASS_lams = np.loadtxt('J_2MASS.TRM')[:,0] * unit['AA']
J_2MASS_trans = np.loadtxt('J_2MASS.TRM')[:,1]
H_2MASS_lams = np.loadtxt('H_2MASS.TRM')[:,0] * unit['AA']
H_2MASS_trans = np.loadtxt('H_2MASS.TRM')[:,1]
K_2MASS_lams = np.loadtxt('K_2MASS.TRM')[:,0] * unit['AA']
K_2MASS_trans = np.loadtxt('K_2MASS.TRM')[:,1]
#Johnson_U_lams = np.loadtxt('johnsonU.datTRM')[:,0] * unit['AA']
#Johnson_U_trans = np.loadtxt('johnsonU.datTRM')[:,1]
Johnson_B_lams = np.loadtxt('johnsonB.datTRM')[:,0] * unit['AA']
Johnson_B_trans = np.loadtxt('johnsonB.datTRM')[:,1]
Johnson_V_lams = np.loadtxt('johnsonV.datTRM')[:,0] * unit['AA']
Johnson_V_trans = np.loadtxt('johnsonV.datTRM')[:,1]
sdss_u_lams = np.loadtxt('sdssu.TRM')[:,0] * unit['AA']
sdss_u_trans = np.loadtxt('sdssu.TRM')[:,1]
sdss_g_lams = np.loadtxt('sdssg.TRM')[:,0] * unit['AA']
sdss_g_trans = np.loadtxt('sdssg.TRM')[:,1]
sdss_r_lams = np.loadtxt('sdssr.TRM')[:,0] * unit['AA']
sdss_r_trans = np.loadtxt('sdssr.TRM')[:,1]
sdss_i_lams = np.loadtxt('sdssi.TRM')[:,0] * unit['AA']
sdss_i_trans = np.loadtxt('sdssi.TRM')[:,1]
sdss_z_lams = np.loadtxt('sdssz.TRM')[:,0] * unit['AA']
sdss_z_trans = np.loadtxt('sdssz.TRM')[:,1]
WISE1_lams = np.loadtxt('Wise1.TRM')[:,0] * unit['AA']
WISE1_trans = np.loadtxt('Wise1.TRM')[:,1]
WISE2_lams = np.loadtxt('Wise2.TRM')[:,0] * unit['AA']
WISE2_trans = np.loadtxt('Wise2.TRM')[:,1]
WISE3_lams = np.loadtxt('Wise3.TRM')[:,0] * unit['AA']
WISE3_trans = np.loadtxt('Wise3.TRM')[:,1]
WISE4_lams = np.loadtxt('Wise4.TRM')[:,0] * unit['AA']
WISE4_trans = np.loadtxt('Wise4.TRM')[:,1]
fuv_galex_lams = np.loadtxt('galex_fuv.TRM')[:,0] * unit['AA']
fuv_galex_trans = np.loadtxt('galex_fuv.TRM')[:,1]
nuv_galex_lams = np.loadtxt('galex_nuv.TRM')[:,0] * unit['AA']
nuv_galex_trans = np.loadtxt('galex_nuv.TRM')[:,1]

#gaia_g_filter = UnitFilter(gaia_g_lams,gaia_g_trans,unit=None)
#J_2MASS_filter = Filter(J_2MASS_lams,J_2MASS_trans,unit='Angstrom')
#H_2MASS_filter = Filter(H_2MASS_lams,H_2MASS_trans,unit='Angstrom')
#K_2MASS_filter = Filter(K_2MASS_lams,K_2MASS_trans,unit='Angstrom')
#Johnson_U_filter = Filter(Johnson_U_lams,Johnson_U_trans,unit='Angstrom')
#Johnson_B_filter = Filter(Johnson_B_lams,Johnson_B_trans,unit='Angstrom',dtype='energy')
#Johnson_V_filter = Filter(Johnson_V_lams,Johnson_V_trans,unit='Angstrom',dtype='energy')
#sdss_u_filter = Filter(sdss_u_lams,sdss_u_trans,unit='Angstrom')
#sdss_g_filter = Filter(sdss_g_lams,sdss_g_trans,unit='Angstrom')
#sdss_r_filter = Filter(sdss_r_lams,sdss_r_trans,unit='Angstrom')
#sdss_i_filter = Filter(sdss_i_lams,sdss_i_trans,unit='Angstrom')
#sdss_z_filter = Filter(sdss_z_lams,sdss_z_trans,unit='Angstrom')
#WISE1_filter = Filter(WISE1_lams,WISE1_trans,unit='Angstrom')
#WISE2_filter = Filter(WISE2_lams,WISE2_trans,unit='Angstrom')
#WISE3_filter = Filter(WISE3_lams,WISE3_trans,unit='Angstrom')
#WISE4_filter = Filter(WISE4_lams,WISE4_trans,unit='Angstrom')
#fuv_galex_filter = Filter(fuv_galex_lams,fuv_galex_trans,unit='Angstrom')
#nuv_galex_filter = Filter(nuv_galex_lams,nuv_galex_trans,unit='Angstrom')

lib = pyphot.get_library()
J_2MASS_filter = lib['2MASS_J']
H_2MASS_filter = lib['2MASS_H']
K_2MASS_filter = lib['2MASS_Ks']
Johnson_B_filter = lib['GROUND_JOHNSON_B']
Johnson_V_filter = lib['GROUND_JOHNSON_V']
sdss_u_filter = lib['SDSS_u']
sdss_g_filter = lib['SDSS_g']
sdss_r_filter = lib['SDSS_r']
sdss_i_filter = lib['SDSS_i']
sdss_z_filter = lib['SDSS_z']
WISE1_filter = lib['WISE_RSR_W1']
WISE2_filter = lib['WISE_RSR_W2']
WISE3_filter = lib['WISE_RSR_W3']
WISE4_filter = lib['WISE_RSR_W4']
fuv_galex_filter = lib['GALEX_FUV']
nuv_galex_filter = lib['GALEX_NUV']

			#read in temp and logg... right now just giving a 10% error budget on temp and logg

os.chdir('/media/intern2/DD82-B92F/RESEARCH_FALL_2022/')
#string_list = []
#with open("Roy_teff_logg_real_names.txt", "r") as file:
#	for line in file:
#		string_value = line.split()[0]
#		string_list.append(string_value)			then just change line.split number and redo for teff, teff_err, logg, and logg_err

#temp = np.loadtxt("Roy_teff_logg_real_names.txt",usecols=0)[26]
#temp_err = np.loadtxt("Roy_teff_logg_real_names.txt",usecols=1)[26]
#logg = np.loadtxt("Roy_teff_logg_real_names.txt",usecols=2)[26]
#logg_err = np.loadtxt("Roy_teff_logg_real_names.txt",usecols=3)[26]


model_directory = '/media/intern2/DD82-B92F/RESEARCH_FALL_2022/stitched_blackbody_binned_models/'
os.chdir(model_directory)
star_name = input("Star Name:")
#temp = np.loadtxt('input_params_for_sed.py')[0]	#gd934 info only
temp = 31041.0
temp_err = 0.01*temp
#logg = np.loadtxt('input_params_for_sed.py')[1]
logg = 5.52
logg_err = 0.01*logg
rv = 3.1
rv_err = 0.1*rv
model_list = os.listdir(model_directory)
numbered_model_list = []
temp_to_save = temp
grav_to_save = logg
temp_err_to_save = temp_err
grav_err_to_save = logg_err
rv_to_save = rv
rv_err_to_save = rv_err


for model in model_list:
	numbered_model_list.append(model[17]+model[18]+model[19]+model[21]+model[22]+model[23])


Temperature_values = []
Temperature_err_values = []
Gravity_values = []
Gravity_err_values = []
EBV_values = []
EBV_err_values = []
Rv_values = []
Rv_err_values = []
Distance_Pc_values = []
Distance_Pc_err_values = []
Radius_values = []
Radius_err_values = []
chi_square_values = []
Mass_values = []
Mass_err_values = []
Absolute_err_values = []
sdss_u_flux_Vals = []
sdss_g_flux_Vals = []
sdss_r_flux_Vals = []
sdss_i_flux_Vals = []
sdss_z_flux_Vals = []
Wise1_flux_Vals = []
Wise2_flux_Vals = []
Wise3_flux_Vals = []
Wise4_flux_Vals = []
Johnson_B_flux_Vals = []
Johnson_V_flux_Vals = []
FUV_flux_Vals = []
NUV_flux_Vals = []
J2MASS_flux_Vals = []
H2MASS_flux_Vals = []
K2MASS_flux_Vals = []
sdss_u_flux_err_Vals = []
sdss_g_flux_err_Vals = []
sdss_r_flux_err_Vals = []
sdss_i_flux_err_Vals = []
sdss_z_flux_err_Vals = []
Wise1_flux_err_Vals = []
Wise2_flux_err_Vals = []
Wise3_flux_err_Vals = []
Wise4_flux_err_Vals = []
Johnson_B_flux_err_Vals = []
Johnson_V_flux_err_Vals = []
FUV_flux_err_Vals = []
NUV_flux_err_Vals = []
J2MASS_flux_err_Vals = []
H2MASS_flux_err_Vals = []
K2MASS_flux_err_Vals = []

G = 6.6743*10**(-11)	#in m^3 kg^-1 s^-2

				#GET COORDINATES FROM STAR NAME
dec_hms = Simbad.query_object(star_name)['DEC']
ra_hms = Simbad.query_object(star_name)['RA']
ra_hms_str = str(ra_hms)
ra_list = [(ra_hms_str[42]+ra_hms_str[43]),(ra_hms_str[45]+ra_hms_str[46]),(ra_hms_str[48]+ra_hms_str[49]+ra_hms_str[50]+ra_hms_str[51]+ra_hms_str[52]+ra_hms_str[53]+ra_hms_str[54])]
hr = float(ra_list[0])
min = float(ra_list[1])
sec = float(ra_list[2])
ra_deg = 15*(hr+(min/60)+(sec/3600))
dec_hms_str = str(dec_hms)
sign = dec_hms_str[42]
dec_list = [(dec_hms_str[42]+dec_hms_str[43]+dec_hms_str[44]),(dec_hms_str[46]+dec_hms_str[47]),(dec_hms_str[49]+dec_hms_str[50]+dec_hms_str[51]+dec_hms_str[52]+dec_hms_str[53]+dec_hms_str[54])]
day = float(dec_list[0])
min = float(dec_list[1])
sec = float(dec_list[2])
if sign == "-":
        dec_deg = day-(min/60)-(sec/3600)
else:
        dec_deg = day+(min/60)+(sec/3600)

co = coords.SkyCoord(ra_deg, dec_deg, unit='deg')
print(co)

                                #E(B-V) ASTROQUERY
dust_radius = 2*u.deg
dust_table = IrsaDust.get_query_table(co, radius=dust_radius, section='ebv')
ebv_temp = dust_table['ext SandF mean']
ebv_err_temp = dust_table['ext SandF std']
ebv = ebv_temp[0]
ebv_err = ebv_err_temp[0]
print('E(B-V) in mag:',ebv)
print('E(B-V) error in mag:',ebv_err)
#print('type E(B-V): ', type(ebv))
#print('type E(B-V) err: ', type(ebv_err))
ebv_to_save = ebv
ebv_err_to_save = ebv_err


				#Rv ASTROQUERY
#Rv_table = dustmaps.(co, radius=dust_radius, section='ebv')
#print(Rv_table[0])

				#SDSS ASTROQUERY
sdss_radius = 0.00005   #in degrees
sdss_query = """ SELECT u,err_u,g,err_g,r,err_r,i,err_i,z,err_z FROM PhotoObj WHERE ra>{ra_min} AND ra<{ra_max} AND dec>{dec_min} AND dec<{dec_max} """.format(ra_min=ra_deg-sdss_radius,ra_max=ra_deg+sdss_radius,dec_min=dec_deg-sdss_radius,dec_max=dec_deg+sdss_radius)
sdss_query_standard = """ SELECT u FROM PhotoObj WHERE ra>18.5775 and ra<18.5776 and dec>-0.8205 and dec<-0.8204 """
sdss_job = SDSS.query_sql(sdss_query,data_release=18)
sdss_job_standard = SDSS.query_sql(sdss_query_standard,data_release=18)

if type(sdss_job) == type(sdss_job_standard):
	u_mag_to_sum,u_mag_err_to_sum,g_mag_to_sum,g_mag_err_to_sum,r_mag_to_sum,r_mag_err_to_sum,i_mag_to_sum,i_mag_err_to_sum,z_mag_to_sum,z_mag_err_to_sum = [],[],[],[],[],[],[],[],[],[]

	u_mag,u_mag_err,g_mag,g_mag_err,r_mag,r_mag_err,i_mag,i_mag_err,z_mag,z_mag_err = sdss_job['u'],sdss_job['err_u'],sdss_job['g'],sdss_job['err_g'],sdss_job['r'],sdss_job['err_r'],sdss_job['i'],sdss_job['err_i'],sdss_job['z'],sdss_job['err_z']

	for num in range(0,len(u_mag)):
        	if u_mag[num]*u_mag_err[num] < u_mag[num]:
                	u_mag_to_sum.append(u_mag[num])
                	u_mag_err_to_sum.append(u_mag_err[num])
        	if g_mag[num]*g_mag_err[num] < g_mag[num]:
                	g_mag_to_sum.append(g_mag[num])
                	g_mag_err_to_sum.append(g_mag_err[num])
        	if r_mag[num]*r_mag_err[num] < r_mag[num]:
                	r_mag_to_sum.append(r_mag[num])
                	r_mag_err_to_sum.append(r_mag_err[num])
        	if i_mag[num]*i_mag_err[num] < i_mag[num]:
                	i_mag_to_sum.append(i_mag[num])
                	i_mag_err_to_sum.append(i_mag_err[num])
        	if z_mag[num]*z_mag_err[num] < z_mag[num]:
                	z_mag_to_sum.append(z_mag[num])
                	z_mag_err_to_sum.append(z_mag_err[num])

	summed_u_mag,summed_u_mag_err,summed_g_mag,summed_g_mag_err,summed_r_mag,summed_r_mag_err,summed_i_mag,summed_i_mag_err,summed_z_mag,summed_z_mag_err = np.sum(u_mag_to_sum),np.sum(u_mag_err_to_sum),np.sum(g_mag_to_sum),np.sum(g_mag_err_to_sum),np.sum(r_mag_to_sum),np.sum(r_mag_err_to_sum),np.sum(i_mag_to_sum),np.sum(i_mag_err_to_sum),np.sum(z_mag_to_sum),np.sum(z_mag_err_to_sum)

	u_mag_avg,u_mag_err_avg,g_mag_avg,g_mag_err_avg,r_mag_avg,r_mag_err_avg,i_mag_avg,i_mag_err_avg,z_mag_avg,z_mag_err_avg = summed_u_mag/(len(u_mag_to_sum)),summed_u_mag_err/(len(u_mag_err_to_sum)),summed_g_mag/(len(g_mag_to_sum)),summed_g_mag_err/(len(g_mag_err_to_sum)),summed_r_mag/(len(r_mag_to_sum)),summed_r_mag_err/(len(r_mag_err_to_sum)),summed_i_mag/(len(i_mag_to_sum)),summed_i_mag_err/(len(i_mag_err_to_sum)),summed_z_mag/(len(z_mag_to_sum)),summed_z_mag_err/(len(z_mag_err_to_sum))

	u_zp_Jy_flux,g_zp_Jy_flux,r_zp_Jy_flux,i_zp_Jy_flux,z_zp_Jy_flux = 3767,3631,3631,3631,3565

	u_flux_avg_Jy,g_flux_avg_Jy,r_flux_avg_Jy,i_flux_avg_Jy,z_flux_avg_Jy = (2*1.4*(1e-10)*u_zp_Jy_flux)*np.sinh((u_mag_avg*np.log(10)/-2.5)-np.log(1.4*(1e-10))),(2*0.9*(1e-10)*g_zp_Jy_flux)*np.sinh((g_mag_avg*np.log(10)/-2.5)-np.log(0.9*(1e-10))),(2*1.2*(1e-10)*r_zp_Jy_flux)*np.sinh((r_mag_avg*np.log(10)/-2.5)-np.log(1.2*(1e-10))),(2*1.8*(1e-10)*i_zp_Jy_flux)*np.sinh((i_mag_avg*np.log(10)/-2.5)-np.log(1.8*(1e-10))),(2*7.4*(1e-10)*z_zp_Jy_flux)*np.sinh((z_mag_avg*np.log(10)/-2.5)-np.log(7.4*(1e-10)))

	err_u_flux_avg_Jy,err_g_flux_avg_Jy,err_r_flux_avg_Jy,err_i_flux_avg_Jy,err_z_flux_avg_Jy = np.abs(((u_flux_avg_Jy*u_mag_err_avg)/(-2.5/np.log(10)))/np.tanh(((u_mag_avg/(-2.5/np.log(10)))-np.log(1.4e-10)))),np.abs(((g_flux_avg_Jy*g_mag_err_avg)/(-2.5/np.log(10)))/np.tanh(((g_mag_avg/(-2.5/np.log(10)))-np.log(0.9e-10)))),np.abs(((r_flux_avg_Jy*r_mag_err_avg)/(-2.5/np.log(10)))/np.tanh(((r_mag_avg/(-2.5/np.log(10)))-np.log(1.2e-10)))),np.abs(((i_flux_avg_Jy*i_mag_err_avg)/(-2.5/np.log(10)))/np.tanh(((i_mag_avg/(-2.5/np.log(10)))-np.log(1.8e-10)))),np.abs(((z_flux_avg_Jy*z_mag_err_avg)/(-2.5/np.log(10)))/np.tanh(((z_mag_avg/(-2.5/np.log(10)))-np.log(7.4e-10))))

	u_lam_eff,g_lam_eff,r_lam_eff,i_lam_eff,z_lam_eff = 3551,4686,6166,7480,8932

	u_flux_avg_erg,g_flux_avg_erg,r_flux_avg_erg,i_flux_avg_erg,z_flux_avg_erg =  (u_flux_avg_Jy*1e-23)/((u_lam_eff**2)/(2.9979e18)),(g_flux_avg_Jy*1e-23)/((g_lam_eff**2)/(2.9979e18)),(r_flux_avg_Jy*1e-23)/((r_lam_eff**2)/(2.9979e18)),(i_flux_avg_Jy*1e-23)/((i_lam_eff**2)/(2.9979e18)),(z_flux_avg_Jy*1e-23)/((z_lam_eff**2)/(2.9979e18))

	err_u_flux_avg_erg,err_g_flux_avg_erg,err_r_flux_avg_erg,err_i_flux_avg_erg,err_z_flux_avg_erg = (err_u_flux_avg_Jy*1e-23)/((u_lam_eff**2)/(2.9979e18)),(err_g_flux_avg_Jy*1e-23)/((g_lam_eff**2)/(2.9979e18)),(err_r_flux_avg_Jy*1e-23)/((r_lam_eff**2)/(2.9979e18)),(err_i_flux_avg_Jy*1e-23)/((i_lam_eff**2)/(2.9979e18)),(err_z_flux_avg_Jy*1e-23)/((z_lam_eff**2)/(2.9979e18))

	print('SDSS:')
	print('u flux in erg', u_flux_avg_erg)
	#print('type sdss u: ', type(u_flux_avg_erg))
	print('g flux in erg', g_flux_avg_erg)
	print('r flux in erg', r_flux_avg_erg)
	print('i flux in erg', i_flux_avg_erg)
	print('z flux in erg', z_flux_avg_erg)
	print('u error flux in erg', err_u_flux_avg_erg)
	#print('type sdss u err: ', type(err_u_flux_avg_erg))
	print('g error flux in erg', err_g_flux_avg_erg)
	print('r error flux in erg', err_r_flux_avg_erg)
	print('i error flux in erg', err_i_flux_avg_erg)
	print('z error flux in erg', err_z_flux_avg_erg)

if type(sdss_job) != type(sdss_job_standard):
	u_flux_avg_erg,g_flux_avg_erg,r_flux_avg_erg,i_flux_avg_erg,z_flux_avg_erg,err_u_flux_avg_erg,err_g_flux_avg_erg,err_r_flux_avg_erg,err_i_flux_avg_erg,err_z_flux_avg_erg = 'NA','NA','NA','NA','NA','NA','NA','NA','NA','NA'
	print('SDSS:')
	print('u flux in erg', u_flux_avg_erg)
	print('g flux in erg', g_flux_avg_erg)
	print('r flux in erg', r_flux_avg_erg)
	print('i flux in erg', i_flux_avg_erg)
	print('z flux in erg', z_flux_avg_erg)
	print('u error flux in erg', err_u_flux_avg_erg)
	print('g error flux in erg', err_g_flux_avg_erg)
	print('r error flux in erg', err_r_flux_avg_erg)
	print('i error flux in erg', err_i_flux_avg_erg)
	print('z error flux in erg', err_z_flux_avg_erg)

sdss_result = SDSS.query_region(co,spectro=True)
sdss_result_type = type(sdss_result)
sdss_str = str(sdss_result_type)
sdss_length = len(sdss_str)

if sdss_length == 35:                                                   #sdss_length=35 corresponds to an astropytable which means sdss_result != None
        sdss_spec = SDSS.get_spectra(matches=sdss_result)[0][1]
        sdss_spec_data = sdss_spec.data
        loglam = sdss_spec_data['loglam']
        SDSS_wavelengths = [10**val for val in loglam]			#observed wavelengths
        SDSS_fluxes = sdss_spec_data['flux']				#observed flux
        #plt.plot(SDSS_wavelengths,SDSS_fluxes)
        #plt.xlabel(r'Wavelength($\AA$)')
        #plt.ylabel(r'Flux(1E-17 erg/cm^2/s/$\AA$)')
        #plt.show()
else:
        print('No Spectrum Available')



                                                #GAIA ASTROQUERY
GAIA_radius = 0.0005    #in degrees
gaia_query = """ SELECT parallax, parallax_error, phot_g_mean_flux, phot_g_mean_flux_error, phot_g_mean_mag,ra,dec FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec),CIRCLE({ra},{dec},{radius}))""".format(ra=ra_deg,dec=dec_deg,radius=GAIA_radius)
gaia_job = Gaia.launch_job(gaia_query)
gaia_results = gaia_job.get_results()
if len(gaia_results) !=0 and gaia_results[0]['parallax'] >= 0:
	print('GAIA:')

	parallax_to_sum,parallax_err_to_sum = [],[]

	parallax,parallax_err = gaia_results['parallax'],gaia_results['parallax_error']

	for val in range(0,len(parallax)):
		if parallax[val]*parallax_err[val] < parallax[val]:
	        	parallax_to_sum.append(parallax[val])
	        	parallax_err_to_sum.append(parallax_err[val])

	summed_parallax,summed_parallax_err = np.sum(parallax_to_sum),np.sum(parallax_err_to_sum)

	avg_parallax,avg_parallax_err = summed_parallax/(len(parallax_to_sum)),summed_parallax_err/(len(parallax_err_to_sum))

	percent_error = avg_parallax_err/avg_parallax*100
	distance = 1/(np.tan(avg_parallax/1000))
	distance_err = percent_error/100*distance
	print('Distance in Pc:',distance)
	#print('distance type: ', type(distance))
	print('Distance error in Pc:',distance_err)
	#print('distance err type: ', type(distance_err))
else:
	distance = 'NA'
	distance_err = 'NA'
	sys.exit("NO DISTANCE MEASUREMENT IN GAIA")		#can change this to a break statement if this is within a for loop of many stars

distance_to_save = distance
distance_err_to_save = distance_err

                                #WISE ASTROQUERY
wise_radius = 0.005*u.deg
wise_job = Irsa.query_region(co, catalog='allwise_p3as_psd', spatial='Cone', radius=wise_radius)
W1_mag,W1_mag_err,W2_mag,W2_mag_err,W3_mag,W3_mag_err,W4_mag,W4_mag_err = wise_job['w1mpro'],wise_job['w1sigmpro'],wise_job['w2mpro'],wise_job['w2sigmpro'],wise_job['w3mpro'],wise_job['w3sigmpro'],wise_job['w4mpro'],wise_job['w4sigmpro']
W1_to_sum,W1_err_to_sum,W2_to_sum,W2_err_to_sum,W3_to_sum,W3_err_to_sum,W4_to_sum,W4_err_to_sum = [],[],[],[],[],[],[],[]

for entry in range(0,len(W1_mag)):
        if type(W1_mag[entry]) == np.float64:
                W1_to_sum.append(W1_mag[entry])

for entry in range(0,len(W1_mag_err)):
        if type(W1_mag_err[entry]) == np.float64:
                W1_err_to_sum.append(W1_mag_err[entry])

for entry in range(0,len(W2_mag)):
        if type(W2_mag[entry]) == np.float64:
                W2_to_sum.append(W2_mag[entry])

for entry in range(0,len(W2_mag_err)):
        if type(W2_mag_err[entry]) == np.float64:
                W2_err_to_sum.append(W2_mag_err[entry])

for entry in range(0,len(W3_mag)):
        if type(W3_mag[entry]) == np.float64:
                W3_to_sum.append(W3_mag[entry])

for entry in range(0,len(W3_mag_err)):
        if type(W3_mag_err[entry]) == np.float64:
                W3_err_to_sum.append(W3_mag_err[entry])

for entry in range(0,len(W4_mag)):
        if type(W4_mag[entry]) == np.float64:
                W4_to_sum.append(W4_mag[entry])

for entry in range(0,len(W4_mag_err)):
        if type(W4_mag_err[entry]) == np.float64:
                W4_err_to_sum.append(W4_mag_err[entry])


if len(W1_to_sum)>=1:
        avg_W1_mag = np.sum(W1_to_sum)/len(W1_to_sum)
else:
        avg_W1_mag = 'NA'

if len(W1_err_to_sum)>=1:
        avg_W1_mag_err = np.sum(W1_err_to_sum)/len(W1_err_to_sum)
else:
        avg_W1_mag_err = 'NA'

if len(W2_to_sum)>=1:
        avg_W2_mag = np.sum(W2_to_sum)/len(W2_to_sum)
else:
        avg_W2_mag = 'NA'
if len(W2_err_to_sum)>=1:
        avg_W2_mag_err = np.sum(W2_err_to_sum)/len(W2_err_to_sum)
else:
        avg_W2_mag_err = 'NA'

if len(W3_to_sum)>=1:
        avg_W3_mag = np.sum(W3_to_sum)/len(W3_to_sum)
else:
        avg_W3_mag = 'NA'

if len(W3_err_to_sum)>=1:
        avg_W3_mag_err = np.sum(W3_err_to_sum)/len(W3_err_to_sum)
else:
        avg_W3_mag_err = 'NA'

if len(W4_to_sum)>=1:
        avg_W4_mag = np.sum(W4_to_sum)/len(W4_to_sum)
else:
        avg_W4_mag = 'NA'

if len(W4_err_to_sum)>=1:
        avg_W4_mag_err = np.sum(W4_err_to_sum)/len(W4_err_to_sum)
else:
        avg_W4_mag_err = 'NA'


W1_zp,W2_zp,W3_zp,W4_zp = 309.540,171.787,31.674,8.363
Jy_list = []
Jy_error_list = []

wise_list = [avg_W1_mag,avg_W1_mag_err,avg_W2_mag,avg_W2_mag_err,avg_W3_mag,avg_W3_mag_err,avg_W4_mag,avg_W4_mag_err]

if type(wise_list[0]) == np.float64:
        if type(wise_list[1]) == np.float64:
                Jy_list.append(W1_zp*10**(-1*avg_W1_mag/2.5))
                Jy_error_list.append(((W1_zp*10**(-1*(avg_W1_mag-avg_W1_mag_err)/2.5))-(W1_zp*10**(-1*(avg_W1_mag+avg_W1_mag_err)/2.5))))
        else:
                Jy_list.append(W1_zp*10**(-1*avg_W1_mag/2.5))
                Jy_error_list.append('NA')

if type(wise_list[0]) != np.float64:
        Jy_list.append('NA')
        Jy_error_list.append('NA')

if type(wise_list[2]) == np.float64:
        if type(wise_list[3]) == np.float64:
                Jy_list.append(W2_zp*10**(-1*avg_W2_mag/2.5))
                Jy_error_list.append(((W2_zp*10**(-1*(avg_W2_mag-avg_W2_mag_err)/2.5))-(W2_zp*10**(-1*(avg_W2_mag+avg_W2_mag_err)/2.5))))
        else:
                Jy_list.append(W2_zp*10**(-1*avg_W2_mag/2.5))
                Jy_error_list.append('NA')

if type(wise_list[2]) != np.float64:
        Jy_list.append('NA')
        Jy_error_list.append('NA')

if type(wise_list[4]) == np.float64:
        if type(wise_list[5]) == np.float64:
                Jy_list.append(W3_zp*10**(-1*avg_W3_mag/2.5))
                Jy_error_list.append(((W3_zp*10**(-1*(avg_W3_mag-avg_W3_mag_err)/2.5))-(W3_zp*10**(-1*(avg_W3_mag+avg_W3_mag_err)/2.5))))
        else:
                Jy_list.append(W3_zp*10**(-1*avg_W3_mag/2.5))
                Jy_error_list.append('NA')

if type(wise_list[4]) != np.float64:
        Jy_list.append('NA')
        Jy_error_list.append('NA')

if type(wise_list[6]) == np.float64:
        if type(wise_list[7]) == np.float64:
                Jy_list.append(W4_zp*10**(-1*avg_W4_mag/2.5))
                Jy_error_list.append(((W4_zp*10**(-1*(avg_W4_mag-avg_W4_mag_err)/2.5))-(W4_zp*10**(-1*(avg_W4_mag+avg_W4_mag_err)/2.5))))
        else:
                Jy_list.append(W4_zp*10**(-1*avg_W4_mag/2.5))
                Jy_error_list.append('NA')

if type(wise_list[6]) != np.float64:
        Jy_list.append('NA')
        Jy_error_list.append('NA')

W1_effwl,W2_effwl,W3_effwl,W4_effwl = (3.368*10**4),(4.618*10**4),(12.082*10**4),(22.194*10**4)
erg_list = []
if type(Jy_list[0]) == np.float64:
        to_erg = Jy_list[0]*10**-23
        band = (W1_effwl**2)/(2.9979*10**18)
        erg_list.append(to_erg/band)
if type(Jy_list[0]) == str:
        erg_list.append('NA')
if type(Jy_list[1]) == np.float64:
        to_erg = Jy_list[1]*10**-23
        band = (W2_effwl**2)/(2.9979*10**18)
        erg_list.append(to_erg/band)
if type(Jy_list[1]) == str:
        erg_list.append('NA')
if type(Jy_list[2]) == np.float64:
        to_erg = Jy_list[2]*10**-23
        band = (W3_effwl**2)/(2.9979*10**18)
        erg_list.append(to_erg/band)
if type(Jy_list[2]) == str:
        erg_list.append('NA')
if type(Jy_list[3]) == np.float64:
        to_erg = Jy_list[3]*10**-23
        band = (W4_effwl**2)/(2.9979*10**18)
        erg_list.append(to_erg/band)
if type(Jy_list[3]) == str:
        erg_list.append('NA')

erg_error_list = []
if type(Jy_error_list[0]) == np.float64:
        to_erg = Jy_error_list[0]*10**-23
        band = (W1_effwl**2)/(2.9979*10**18)
        erg_error_list.append(to_erg/band)
if type(Jy_error_list[0]) == str:
        erg_error_list.append('NA')
if type(Jy_error_list[1]) == np.float64:
        to_erg = Jy_error_list[1]*10**-23
        band = (W2_effwl**2)/(2.9979*10**18)
        erg_error_list.append(to_erg/band)
if type(Jy_error_list[1]) == str:
        erg_error_list.append('NA')
if type(Jy_error_list[2]) == np.float64:
        to_erg = Jy_error_list[2]*10**-23
        band = (W3_effwl**2)/(2.9979*10**18)
        erg_error_list.append(to_erg/band)
if type(Jy_error_list[2]) == str:
        erg_error_list.append('NA')
if type(Jy_error_list[3]) == np.float64:
        to_erg = Jy_error_list[3]*10**-23
        band = (W4_effwl**2)/(2.9979*10**18)
        erg_error_list.append(to_erg/band)
if type(Jy_error_list[3]) == str:
        erg_error_list.append('NA')

print('WISE:')
print('W1 erg:',erg_list[0])
#print('type W1: ',type(erg_list[0]))
print('W1 error erg:', erg_error_list[0])
#print('type W1 err: ', type(erg_error_list[0]))
print('W2 erg:',erg_list[1])
print('W2 error erg:', erg_error_list[1])
print('W3 erg:',erg_list[2])
print('W3 error erg:', erg_error_list[2])
print('W4 erg:',erg_list[3])
print('W4 error erg:', erg_error_list[3])

                        #GALEX ASTROQUERY
galex_region = Catalogs.query_region(co,radius=0.002,catalog="galex")
galex_ra = galex_region['ra']
galex_dec = galex_region['dec']
galex_fuv_mag = galex_region['fuv_mag']
galex_fuv_mag_err = galex_region['fuv_magerr']
galex_nuv_mag = galex_region['nuv_mag']
galex_nuv_mag_err = galex_region['nuv_magerr']
galex_fuv_flux = galex_region['fuv_flux']
galex_fuv_flux_err = galex_region['fuv_fluxerr']
galex_nuv_flux = galex_region['nuv_flux']
galex_nuv_flux_err = galex_region['nuv_fluxerr']

if len(galex_region)!=0:
	sum_ra,sum_dec,sum_fuv_mag,sum_fuv_mag_err,sum_fuv_flux,sum_fuv_flux_err,sum_nuv_mag,sum_nuv_mag_err,sum_nuv_flux,sum_nuv_flux_err = [],[],[],[],[],[],[],[],[],[]

	for n in range(0,len(galex_ra)):
		if type(galex_ra[n]) == np.float64:
                	sum_ra.append(galex_ra[n])
		if type(galex_ra[n]) != np.float64:
			continue
		if type(galex_dec[n]) == np.float64:
			sum_dec.append(galex_dec[n])
		if type(galex_dec[n]) != np.float64:
			continue
		if type(galex_fuv_mag[n]) == np.float64:
			sum_fuv_mag.append(galex_fuv_mag[n])
		if type(galex_fuv_mag[n]) != np.float64:
			continue
		if type(galex_fuv_mag_err[n]) == np.float64:
			sum_fuv_mag_err.append(galex_fuv_mag_err[n])
		if type(galex_fuv_mag_err[n]) != np.float64:
			continue
		if type(galex_fuv_flux[n]) == np.float64:
			sum_fuv_flux.append(galex_fuv_flux[n])
		if type(galex_fuv_flux[n]) != np.float64:
			continue
		if type(galex_fuv_flux_err[n]) == np.float64:
			sum_fuv_flux_err.append(galex_fuv_flux_err[n])
		if type(galex_fuv_flux_err[n]) != np.float64:
			continue
		if type(galex_nuv_mag[n]) == np.float64:
			sum_nuv_mag.append(galex_nuv_mag[n])
		if type(galex_nuv_mag[n]) != np.float64:
			continue
		if type(galex_nuv_mag_err[n]) == np.float64:
			sum_nuv_mag_err.append(galex_nuv_mag_err[n])
		if type(galex_nuv_mag_err[n]) != np.float64:
			continue
		if type(galex_nuv_flux[n]) == np.float64:
			sum_nuv_flux.append(galex_nuv_flux[n])
		if type(galex_nuv_flux[n]) != np.float64:
			continue
		if type(galex_nuv_flux_err[n]) == np.float64:
			sum_nuv_flux_err.append(galex_nuv_flux_err[n])
		if type(galex_nuv_flux_err[n]) != np.float64:
			continue

	if len(sum_ra)!=0:
		avg_galex_ra = np.sum(sum_ra)/len(sum_ra)
	if len(sum_ra)==0:
		avg_galex_ra = 'NA'
	if len(sum_dec)!=0:
		avg_galex_dec = np.sum(sum_dec)/len(sum_dec)
	if len(sum_dec)==0:
		avg_galex_dec = 'NA'
	if len(sum_fuv_mag)!=0:
		avg_galex_fuv_mag = np.sum(sum_fuv_mag)/len(sum_fuv_mag)
	if len(sum_fuv_mag)==0:
		avg_galex_fuv_mag = 'NA'
	if len(sum_fuv_mag_err)!=0:
		avg_galex_fuv_mag_err = np.sum(sum_fuv_mag_err)/len(sum_fuv_mag_err)
	if len(sum_fuv_mag_err)==0:
		avg_galex_fuv_mag_err = 'NA'
	if len(sum_fuv_flux)!=0:
		avg_galex_fuv_flux = np.sum(sum_fuv_flux)/len(sum_fuv_flux)
	if len(sum_fuv_flux)==0:
		avg_galex_fuv_flux = 'NA'
	if len(sum_fuv_flux_err)!=0:
		avg_galex_fuv_flux_err = np.sum(sum_fuv_flux_err)/len(sum_fuv_flux_err)
	if len(sum_fuv_flux_err)==0:
		avg_galex_fuv_flux_err = 'NA'
	if len(sum_nuv_mag)!=0:
		avg_galex_nuv_mag = np.sum(sum_nuv_mag)/len(sum_nuv_mag)
	if len(sum_nuv_mag)==0:
		avg_galex_nuv_mag = 'NA'
	if len(sum_nuv_mag_err)!=0:
		avg_galex_nuv_mag_err = np.sum(sum_nuv_mag_err)/len(sum_nuv_mag_err)
	if len(sum_nuv_mag_err)==0:
		avg_galex_nuv_mag_err = 'NA'
	if len(sum_nuv_flux)!=0:
		avg_galex_nuv_flux = np.sum(sum_nuv_flux)/len(sum_nuv_flux)
	if len(sum_nuv_flux)==0:
		avg_galex_nuv_flux = 'NA'
	if len(sum_nuv_flux_err)!=0:
		avg_galex_nuv_flux_err = np.sum(sum_nuv_flux_err)/len(sum_nuv_flux_err)
	if len(sum_nuv_flux_err)==0:
		avg_galex_nuv_flux_err = 'NA'


	if type(avg_galex_fuv_flux) != str:
		fuv_micro_jansky = avg_galex_fuv_flux
		fuv_jansky = fuv_micro_jansky/(10**6)
		fuv_temp = fuv_jansky*(10**(-23))
		fuv_erg = fuv_temp/((1516**2)/(2.9979*10**18))
	if type(avg_galex_fuv_flux) == str:
		fuv_erg = 'NA'
	if type(avg_galex_fuv_flux_err) != str:
		fuv_err_micro_jansky = avg_galex_fuv_flux + avg_galex_fuv_flux_err		#fuv err added to fuv flux
		fuv_err_jansky = fuv_err_micro_jansky/(10**6)
		fuv_err_temp = fuv_err_jansky*(10**(-23))
		fuv_err_temp2 = fuv_err_temp/((1516**2)/(2.9979*10**18))
		fuv_err_erg = fuv_err_temp2-fuv_erg
	if type(avg_galex_fuv_flux_err) == str:
		fuv_err_erg = 'NA'
	if type(avg_galex_nuv_flux) != str:
		nuv_micro_jansky = avg_galex_nuv_flux
		nuv_jansky = nuv_micro_jansky/(10**6)
		nuv_temp = nuv_jansky*(10**(-23))
		nuv_erg = nuv_temp/((2267**2)/(2.9979*10**18))
	if type(avg_galex_nuv_flux) == str:
		nuv_erg = 'NA'
	if type(avg_galex_nuv_flux_err) != str:
		nuv_err_micro_jansky = avg_galex_nuv_flux + avg_galex_nuv_flux_err
		nuv_err_jansky = nuv_err_micro_jansky/(10**6)
		nuv_err_temp = nuv_err_jansky*(10**(-23))
		nuv_err_temp2 = nuv_err_temp/((2267**2)/(2.9979*10**18))
		nuv_err_erg = nuv_err_temp2-nuv_erg
	if type(avg_galex_nuv_flux_err) == str:
		nuv_err_erg = 'NA'

else:
	fuv_erg = 'NA'
	fuv_err_erg = 'NA'
	nuv_erg = 'NA'
	nuv_err_erg = 'NA'

print('GALEX:')
print('fuv erg:',fuv_erg)
#print('type fuv: ', type(fuv_erg))
print('fuv error erg:',fuv_err_erg)
#print('type fuv err: ', type(fuv_err_erg))
print('nuv erg:',nuv_erg)
print('nuv error erg:',nuv_err_erg)


						#2MASS ASTROQUERY
vizier_radius = 0.01*u.deg

vizier_result = Vizier.query_region(co,radius=vizier_radius,catalog='II/246')
vizier_result_str = str(Vizier.query_region(co,radius=vizier_radius,catalog='II/246'))
if vizier_result_str[0] != 'E':
	Jmag_list = vizier_result[0]['Jmag']
	Jmag_err_list = vizier_result[0]['e_Jmag']
	Hmag_list = vizier_result[0]['Hmag']
	Hmag_err_list = vizier_result[0]['e_Hmag']
	Kmag_list = vizier_result[0]['Kmag']
	Kmag_err_list = vizier_result[0]['e_Kmag']

	float64_Jmag_list = []
	float64_Jmag_err_list = []
	float64_Hmag_list = []
	float64_Hmag_err_list = []
	float64_Kmag_list = []
	float64_Kmag_err_list = []

	for entry in Jmag_list:
	        float64_Jmag_list.append(np.float64(entry))
	for entry in Hmag_list:
	        float64_Hmag_list.append(np.float64(entry))
	for entry in Kmag_list:
	        float64_Kmag_list.append(np.float64(entry))

	for entry in Jmag_err_list:
	        if type(entry) == np.float32:
	                float64_Jmag_err_list.append(np.float64(entry))
	        else:
	                float64_Jmag_err_list.append('NA')
	for entry in Hmag_err_list:
	        if type(entry) == np.float32:
	                float64_Hmag_err_list.append(np.float64(entry))
	        else:
	                float64_Hmag_err_list.append('NA')
	for entry in Kmag_err_list:
	        if type(entry) == np.float32:
	                float64_Kmag_err_list.append(np.float64(entry))
	        else:
	                float64_Kmag_err_list.append('NA')


	Jmags_to_sum = []
	Jmag_errs_to_sum = []
	Hmags_to_sum = []
	Hmag_errs_to_sum = []
	Kmags_to_sum = []
	Kmag_errs_to_sum = []
	for i in range(0,len(float64_Jmag_list)):
	        if type(float64_Jmag_err_list[i]) != str:
	                Jmags_to_sum.append(float64_Jmag_list[i])
	                Jmag_errs_to_sum.append(float64_Jmag_err_list[i])
	for i in range(0,len(float64_Hmag_list)):
	        if type(float64_Hmag_err_list[i]) != str:
	                Hmags_to_sum.append(float64_Hmag_list[i])
	                Hmag_errs_to_sum.append(float64_Hmag_err_list[i])
	for i in range(0,len(float64_Kmag_list)):
	        if type(float64_Kmag_err_list[i]) != str:
	                Kmags_to_sum.append(float64_Kmag_list[i])
	                Kmag_errs_to_sum.append(float64_Kmag_err_list[i])

	if len(Jmags_to_sum)!=0:
		avg_Jmag = np.sum(Jmags_to_sum)/len(Jmags_to_sum)
	if len(Jmags_to_sum) == 0:
		avg_Jmag = 'NA'
	if len(Jmag_errs_to_sum)!=0:
		avg_Jmag_err = np.sum(Jmag_errs_to_sum)/len(Jmag_errs_to_sum)
	if len(Jmag_errs_to_sum)==0:
		avg_Jmag_err = 'NA'
	if len(Hmags_to_sum)!=0:
		avg_Hmag = np.sum(Hmags_to_sum)/len(Hmags_to_sum)
	if len(Hmags_to_sum)==0:
		avg_Hmag = 'NA'
	if len(Hmag_errs_to_sum)!=0:
		avg_Hmag_err = np.sum(Hmag_errs_to_sum)/len(Hmag_errs_to_sum)
	if len(Hmag_errs_to_sum)==0:
		avg_Hmag_err = 'NA'
	if len(Kmags_to_sum)!=0:
		avg_Kmag = np.sum(Kmags_to_sum)/len(Kmags_to_sum)
	if len(Kmags_to_sum)==0:
		avg_Kmag = 'NA'
	if len(Kmag_errs_to_sum)!=0:
		avg_Kmag_err = np.sum(Kmag_errs_to_sum)/len(Kmag_errs_to_sum)
	if len(Kmag_errs_to_sum)==0:
		avg_Kmag_err = 'NA'


	if type(avg_Jmag)!=str:
		vizier_J_flux_erg = 3.129*(10**-10)*(10**(avg_Jmag/-2.512))
	else:
		vizier_J_flux_erg = 'NA'
	if type(avg_Hmag)!=str:
		vizier_H_flux_erg = 1.133*(10**-10)*(10**(avg_Hmag/-2.512))
	else:
		vizier_H_flux_erg = 'NA'
	if type(avg_Kmag)!=str:
		vizier_K_flux_erg = 4.283*(10**-11)*(10**(avg_Kmag/-2.512))
	else:
		vizier_K_flux_erg = 'NA'



	if type(avg_Jmag)!=str and type(avg_Jmag_err)!=str:
		to_get_Jerr = avg_Jmag-avg_Jmag_err
		temp_Jerr_flux = 3.129*(10**-10)*(10**(to_get_Jerr/-2.512))
		vizier_J_flux_err_erg = temp_Jerr_flux-vizier_J_flux_erg
	else:
		vizier_J_flux_err_erg = 'NA'
	if type(avg_Hmag)!=str and type(avg_Hmag_err)!=str:
		to_get_Herr = avg_Hmag-avg_Hmag_err
		temp_Herr_flux = 1.133*(10**-10)*(10**(to_get_Herr/-2.512))
		vizier_H_flux_err_erg = temp_Herr_flux-vizier_H_flux_erg
	else:
		vizier_H_flux_err_erg = 'NA'
	if type(avg_Kmag)!=str and type(avg_Kmag_err)!=str:
		to_get_Kerr = avg_Kmag-avg_Kmag_err
		temp_Kerr_flux = 4.283*(10**-11)*(10**(to_get_Kerr/-2.512))
		vizier_K_flux_err_erg = temp_Kerr_flux-vizier_K_flux_erg
	else:
		vizier_K_flux_err_erg = 'NA'

	print('2MASS:')
	print('2MASS J flux erg:',vizier_J_flux_erg)
	#print('type J 2MASS: ', type(vizier_J_flux_erg)) 
	print('2MASS J flux error erg:',vizier_J_flux_err_erg)
	#print('type J 2MASS err: ', type(vizier_J_flux_err_erg))
	print('2MASS H flux erg:',vizier_H_flux_erg)
	print('2MASS H flux error erg:',vizier_H_flux_err_erg)
	print('2MASS K flux erg:',vizier_K_flux_erg)
	print('2MASS K flux error erg:',vizier_K_flux_err_erg)

else:
	vizier_J_flux_erg = 'NA'
	vizier_J_flux_err_erg = 'NA'
	vizier_H_flux_erg = 'NA'
	vizier_H_flux_err_erg = 'NA'
	vizier_K_flux_erg = 'NA'
	vizier_K_flux_err_erg = 'NA'
	print('2MASS:')
	print('2MASS J flux erg:','NA')
	print('2MASS J flux error erg:','NA')
	print('2MASS H flux erg:','NA')
	print('2MASS H flux error erg:','NA')
	print('2MASS K flux erg:','NA')
	print('2MASS K flux error erg:','NA')


			#Johnson B and V Astroquery
if  err_u_flux_avg_erg == 'NA' or err_g_flux_avg_erg == 'NA':
	Vizier.ROW_LIMIT = -1
	Vizier.columns = ['all']
	Johnson_result = Vizier.query_object(star_name,catalog='J/A+A/600/A50')
	Johnson_result_str = str(Vizier.query_object(star_name,catalog='J/A+A/600/A50'))
	if Johnson_result_str[0] != 'E':
		check_if_data_present = Johnson_result[0]['BAPASS']
		if type(check_if_data_present[0]) == np.float32:
			Johnson_B_mag = np.float64(Johnson_result[0]['BAPASS'])
			Johnson_B_err_mag = np.float64(Johnson_result[0]['e_BAPASS'])
			Johnson_V_mag = np.float64(Johnson_result[0]['VAPASS'])
			Johnson_V_err_mag = np.float64(Johnson_result[0]['e_VAPASS'])
			Johnson_B_flux_erg = 6.32*(10**-9)*(10**(Johnson_B_mag/-2.512))
			Johnson_V_flux_erg = 3.631*(10**-9)*(10**(Johnson_V_mag/-2.512))

			to_get_B_err = Johnson_B_mag-Johnson_B_err_mag
			temp_B_flux = 6.32*(10**-9)*(10**(to_get_B_err/-2.512))
			Johnson_B_flux_err_erg = temp_B_flux-Johnson_B_flux_erg
			to_get_V_err = Johnson_V_mag-Johnson_V_err_mag
			temp_V_flux = 3.631*(10**-9)*(10**(to_get_V_err/-2.512))
			Johnson_V_flux_err_erg = temp_V_flux-Johnson_V_flux_erg
			print('JOHNSON:')
			print('Johnson B flux erg: ',Johnson_B_flux_erg)
			#print('type Johnson B: ', type(Johnson_B_flux_erg))
			print('Johnson B flux err erg: ',Johnson_B_flux_err_erg)
			#print('type Johnson B err: ', type(Johnson_B_flux_err_erg))
			print('Johnson V flux erg: ',Johnson_V_flux_erg)
			print('Johnson V flux err erg: ',Johnson_V_flux_err_erg)
		else:
			Johnson_B_flux_erg = 'NA'
			Johnson_B_flux_err_erg = 'NA'
			Johnson_V_flux_erg = 'NA'
			Johnson_V_flux_err_erg = 'NA'
			print('JOHNSON:')
			print('NO APASS DATA IN GEIERS DATABASE')

	else:
		Johnson_B_flux_erg = 'NA'
		Johnson_B_flux_err_erg = 'NA'
		Johnson_V_flux_erg = 'NA'
		Johnson_V_flux_err_erg = 'NA'
		print('JOHNSON:')
		print('THIS STAR DOES NOT EXIST IN GEIERS DATABASE')

if  err_u_flux_avg_erg != 'NA' or err_g_flux_avg_erg != 'NA':
	Johnson_B_flux_erg = 'NA'
	Johnson_B_flux_err_erg = 'NA'
	Johnson_V_flux_erg = 'NA'
	Johnson_V_flux_err_erg = 'NA'
	print('JOHNSON:')
	print('ALREADY HAVE SDSS DATA')







							#START OF MODEL FITTING
inputs = [star_name,temp,temp_err,logg,logg_err,ebv,ebv_err,u_flux_avg_erg,err_u_flux_avg_erg,r_flux_avg_erg,err_r_flux_avg_erg,g_flux_avg_erg,err_g_flux_avg_erg,i_flux_avg_erg,err_i_flux_avg_erg,z_flux_avg_erg,err_z_flux_avg_erg,distance,distance_err,erg_list[0],erg_error_list[0],erg_list[1],erg_error_list[1],erg_list[2],erg_error_list[2],erg_list[3],erg_error_list[3],fuv_erg,fuv_err_erg,nuv_erg,nuv_err_erg,rv,rv_err,vizier_J_flux_erg,vizier_J_flux_err_erg,vizier_H_flux_erg,vizier_H_flux_err_erg,vizier_K_flux_erg,vizier_K_flux_err_erg,Johnson_B_flux_erg,Johnson_B_flux_err_erg,Johnson_V_flux_erg,Johnson_V_flux_err_erg]

for i in range(0,10):				#number of models to be generated
	rand_for_plus_minus = random.rand(1)
	rand_err_for_Teff = (random.rand(1))*inputs[2]
	rand_err_for_LOGg = (random.rand(1))*inputs[4]
	rand_err_for_EBV = (random.rand(1))*inputs[6]
	rand_err_for_Rv = (random.rand(1))*inputs[32]
	if type(inputs[8]) == np.float64:
		rand_err_for_sdss_u_flux = (random.rand(1))*inputs[8]
	else:
		rand_err_for_sdss_u_flux = 'NA'
	if type(inputs[10]) == np.float64:
		rand_err_for_sdss_g_flux = (random.rand(1))*inputs[10]
	else:
		rand_err_for_sdss_g_flux = 'NA'
	if type(inputs[12]) == np.float64:
		rand_err_for_sdss_r_flux = (random.rand(1))*inputs[12]
	else:
		rand_err_for_sdss_r_flux = 'NA'
	if type(inputs[14]) == np.float64:
		rand_err_for_sdss_i_flux = (random.rand(1))*inputs[14]
	else:
		rand_err_for_sdss_i_flux = 'NA'
	if type(inputs[16]) == np.float64:
		rand_err_for_sdss_z_flux = (random.rand(1))*inputs[16]
	else:
		rand_err_for_sdss_z_flux = 'NA'
	if type(inputs[18]) == np.float64:
		rand_err_for_DISTANCE = (random.rand(1))*inputs[18]
	else:
		rand_err_for_DISTANCE = 'NA'
	if type(inputs[20]) == np.float64:
		rand_err_for_W1_flux = (random.rand(1))*inputs[20]
	else:
		rand_err_for_W1_flux = 'NA'
	if type(inputs[22]) == np.float64:
		rand_err_for_W2_flux = (random.rand(1))*inputs[22]
	else:
		rand_err_for_W2_flux = 'NA'
	if type(inputs[24]) == np.float64:
		rand_err_for_W3_flux = (random.rand(1))*inputs[24]
	else:
		rand_err_for_W3_flux = 'NA'
	if type(inputs[26]) == np.float64:
		rand_err_for_W4_flux = (random.rand(1))*inputs[26]
	else:
		rand_err_for_W4_flux = 'NA'

	if type(inputs[28]) == np.float64:
		rand_err_for_FUV_flux = (random.rand(1))*inputs[28]
	else:
		rand_err_for_FUV_flux = 'NA'
	if type(inputs[30]) == np.float64:
		rand_err_for_NUV_flux = (random.rand(1))*inputs[30]
	else:
		rand_err_for_NUV_flux = 'NA'

	if type(inputs[34]) == np.float64:
		rand_err_for_J_flux_2MASS = (random.rand(1))*inputs[34]
	else:
		rand_err_for_J_flux_2MASS = 'NA'

	if type(inputs[36]) == np.float64:
		rand_err_for_H_flux_2MASS = (random.rand(1))*inputs[36]
	else:
		rand_err_for_H_flux_2MASS = 'NA'

	if type(inputs[38]) == np.float64:
		rand_err_for_K_flux_2MASS = (random.rand(1))*inputs[38]
	else:
		rand_err_for_K_flux_2MASS = 'NA'
	if type(inputs[40]) == np.float64:
		rand_err_for_Johnson_B_flux = (random.rand(1))*inputs[40]
	else:
		rand_err_for_Johnson_B_flux = 'NA'
	if type(inputs[42]) == np.float64:
		rand_err_for_Johnson_V_flux = (random.rand(1))*inputs[42]
	else:
		rand_err_for_Johnson_V_flux = 'NA'


	if rand_for_plus_minus <= 0.5:
		Teff = inputs[1]+rand_err_for_Teff
		Teff = float(Teff)
		LOGg = inputs[3]+rand_err_for_LOGg
		LOGg = float(LOGg)
		EBV = inputs[5]+rand_err_for_EBV
		EBV = float(EBV)
		Rv = inputs[31]+rand_err_for_Rv
		Rv = float(Rv)
		if type(rand_err_for_sdss_u_flux) != str:
			sdss_u_flux = inputs[7]+rand_err_for_sdss_u_flux
			sdss_u_flux = float(sdss_u_flux)
		else:
			sdss_u_flux = 'NA'
		if type(rand_err_for_sdss_g_flux) != str:
			sdss_g_flux = inputs[9]+rand_err_for_sdss_g_flux
			sdss_g_flux = float(sdss_g_flux)
		else:
			sdss_g_flux = 'NA'
		if type(rand_err_for_sdss_r_flux) != str:
			sdss_r_flux = inputs[11]+rand_err_for_sdss_r_flux
			sdss_r_flux = float(sdss_r_flux)
		else:
			sdss_r_flux = 'NA'
		if type(rand_err_for_sdss_i_flux) != str:
			sdss_i_flux = inputs[13]+rand_err_for_sdss_i_flux
			sdss_i_flux = float(sdss_i_flux)
		else:
			sdss_i_flux = 'NA'
		if type(rand_err_for_sdss_z_flux) != str:
			sdss_z_flux = inputs[15]+rand_err_for_sdss_z_flux
			sdss_z_flux = float(sdss_z_flux)
		else:
			sdss_z_flux = 'NA'

		if type(rand_err_for_DISTANCE) != str:
			DISTANCE = inputs[17]+rand_err_for_DISTANCE
			DISTANCE = float(DISTANCE)
		else:
			DISTANCE = 'NA'

		if type(rand_err_for_W1_flux) != str:
			W1_flux = inputs[19]+rand_err_for_W1_flux
			W1_flux = float(W1_flux)
		else:
			W1_flux = 'NA'
		if type(rand_err_for_W2_flux) != str:
			W2_flux = inputs[21]+rand_err_for_W2_flux
			W2_flux = float(W2_flux)
		else:
			W2_flux = 'NA'
		if type(rand_err_for_W3_flux) != str:
			W3_flux = inputs[23]+rand_err_for_W3_flux
			W3_flux = float(W3_flux)
		else:
			W3_flux = 'NA'
		if type(rand_err_for_W4_flux) != str:
			W4_flux = inputs[25]+rand_err_for_W4_flux
			W4_flux = float(W4_flux)
		else:
			W4_flux = 'NA'
		if type(rand_err_for_FUV_flux) != str:
			FUV_flux = inputs[27]+rand_err_for_FUV_flux
			FUV_flux = float(FUV_flux)
		else:
			FUV_flux = 'NA'
		if type(rand_err_for_NUV_flux) != str:
			NUV_flux = inputs[29]+rand_err_for_NUV_flux
			NUV_flux = float(NUV_flux)
		else:
			NUV_flux = 'NA'
		if type(rand_err_for_J_flux_2MASS) != str:
			J_flux_2MASS = inputs[33]+rand_err_for_J_flux_2MASS
			J_flux_2MASS = float(J_flux_2MASS)
		else:
			J_flux_2MASS = 'NA'
		if type(rand_err_for_H_flux_2MASS) != str:
			H_flux_2MASS = inputs[35]+rand_err_for_H_flux_2MASS
			H_flux_2MASS = float(H_flux_2MASS)
		else:
			H_flux_2MASS = 'NA'
		if type(rand_err_for_K_flux_2MASS) != str:
			K_flux_2MASS = inputs[37]+rand_err_for_K_flux_2MASS
			K_flux_2MASS = float(K_flux_2MASS)
		else:
			K_flux_2MASS = 'NA'
		if type(rand_err_for_Johnson_B_flux) != str:
			Johnson_B_flux = inputs[39]+rand_err_for_Johnson_B_flux
			Johnson_B_flux = float(Johnson_B_flux)
		else:
			Johnson_B_flux = 'NA'
		if type(rand_err_for_Johnson_V_flux) != str:
			Johnson_V_flux = inputs[41]+rand_err_for_Johnson_V_flux
			Johnson_V_flux = float(Johnson_V_flux)
		else:
			Johnson_V_flux = 'NA'


	else:
		Teff = inputs[1]-rand_err_for_Teff
		Teff = float(Teff)
		LOGg = inputs[3]-rand_err_for_LOGg
		LOGg = float(LOGg)
		EBV = inputs[5]-rand_err_for_EBV
		EBV = float(EBV)
		Rv = inputs[31]-rand_err_for_Rv
		Rv = float(Rv)
		if type(rand_err_for_sdss_u_flux) != str:
			sdss_u_flux = inputs[7]-rand_err_for_sdss_u_flux
			sdss_u_flux = float(sdss_u_flux)
		else:
			sdss_u_flux = 'NA'
		if type(rand_err_for_sdss_g_flux) != str:
			sdss_g_flux = inputs[9]-rand_err_for_sdss_g_flux
			sdss_g_flux = float(sdss_g_flux)
		else:
			sdss_g_flux = 'NA'
		if type(rand_err_for_sdss_r_flux) != str:
			sdss_r_flux = inputs[11]-rand_err_for_sdss_r_flux
			sdss_r_flux = float(sdss_r_flux)
		else:
			sdss_r_flux = 'NA'
		if type(rand_err_for_sdss_i_flux) != str:
			sdss_i_flux = inputs[13]-rand_err_for_sdss_i_flux
			sdss_i_flux = float(sdss_i_flux)
		else:
			sdss_i_flux = 'NA'
		if type(rand_err_for_sdss_z_flux) != str:
			sdss_z_flux = inputs[15]-rand_err_for_sdss_z_flux
			sdss_z_flux = float(sdss_z_flux)
		else:
			sdss_z_flux = 'NA'

		if type(rand_err_for_DISTANCE) != str:
			DISTANCE = inputs[17]-rand_err_for_DISTANCE
			DISTANCE = float(DISTANCE)
		else:
			DISTANCE = 'NA'

		if type(rand_err_for_W1_flux) != str:
			W1_flux = inputs[19]-rand_err_for_W1_flux
			W1_flux = float(W1_flux)
		else:
			W1_flux = 'NA'
		if type(rand_err_for_W2_flux) != str:
			W2_flux = inputs[21]-rand_err_for_W2_flux
			W2_flux = float(W2_flux)
		else:
			W2_flux = 'NA'
		if type(rand_err_for_W3_flux) != str:
			W3_flux = inputs[23]-rand_err_for_W3_flux
			W3_flux = float(W3_flux)
		else:
			W3_flux = 'NA'
		if type(rand_err_for_W4_flux) != str:
			W4_flux = inputs[25]-rand_err_for_W4_flux
			W4_flux = float(W4_flux)
		else:
			W4_flux = 'NA'
		if type(rand_err_for_FUV_flux) != str:
			FUV_flux = inputs[27]-rand_err_for_FUV_flux
			FUV_flux = float(FUV_flux)
		else:
			FUV_flux = 'NA'
		if type(rand_err_for_NUV_flux) != str:
			NUV_flux = inputs[29]-rand_err_for_NUV_flux
			NUV_flux = float(NUV_flux)
		else:
			NUV_flux = 'NA'
		if type(rand_err_for_J_flux_2MASS) != str:
			J_flux_2MASS = inputs[33]-rand_err_for_J_flux_2MASS
			J_flux_2MASS = float(J_flux_2MASS)
		else:
			J_flux_2MASS = 'NA'
		if type(rand_err_for_H_flux_2MASS) != str:
			H_flux_2MASS = inputs[35]-rand_err_for_H_flux_2MASS
			H_flux_2MASS = float(H_flux_2MASS)
		else:
			H_flux_2MASS = 'NA'
		if type(rand_err_for_K_flux_2MASS) != str:
			K_flux_2MASS = inputs[37]-rand_err_for_K_flux_2MASS
			K_flux_2MASS = float(K_flux_2MASS)
		else:
			K_flux_2MASS = 'NA'
		if type(rand_err_for_Johnson_B_flux) != str:
			Johnson_B_flux = inputs[39]-rand_err_for_Johnson_B_flux
			Johnson_B_flux = float(Johnson_B_flux)
		else:
			Johnson_B_flux = 'NA'
		if type(rand_err_for_Johnson_V_flux) != str:
			Johnson_V_flux = inputs[41]-rand_err_for_Johnson_V_flux
			Johnson_V_flux = float(Johnson_V_flux)
		else:
			Johnson_V_flux = 'NA'


	Temperature_values.append(Teff)
	#Temperature_err_values.append(rand_err_for_Teff)
	Gravity_values.append(LOGg)
	#Gravity_err_values.append(rand_err_for_LOGg)
	EBV_values.append(EBV)
	#EBV_err_values.append(rand_err_for_EBV)
	Rv_values.append(Rv)
	#Rv_err_values.append(rand_err_for_Rv)
	Distance_Pc_values.append(DISTANCE)
	#Distance_Pc_err_values.append(rand_err_for_DISTANCE)
	sdss_u_flux_Vals.append(sdss_u_flux)
	sdss_g_flux_Vals.append(sdss_g_flux)
	sdss_r_flux_Vals.append(sdss_r_flux)
	sdss_i_flux_Vals.append(sdss_i_flux)
	sdss_z_flux_Vals.append(sdss_z_flux)
	Wise1_flux_Vals.append(W1_flux)
	Wise2_flux_Vals.append(W2_flux)
	Wise3_flux_Vals.append(W3_flux)
	Wise4_flux_Vals.append(W4_flux)
	Johnson_B_flux_Vals.append(Johnson_B_flux)
	Johnson_V_flux_Vals.append(Johnson_V_flux)
	FUV_flux_Vals.append(FUV_flux)
	NUV_flux_Vals.append(NUV_flux)
	J2MASS_flux_Vals.append(J_flux_2MASS)
	H2MASS_flux_Vals.append(H_flux_2MASS)
	K2MASS_flux_Vals.append(K_flux_2MASS)
	sdss_u_flux_err_Vals.append(rand_err_for_sdss_u_flux)
	sdss_g_flux_err_Vals.append(rand_err_for_sdss_g_flux)
	sdss_r_flux_err_Vals.append(rand_err_for_sdss_r_flux)
	sdss_i_flux_err_Vals.append(rand_err_for_sdss_i_flux)
	sdss_z_flux_err_Vals.append(rand_err_for_sdss_z_flux)
	Wise1_flux_err_Vals.append(rand_err_for_W1_flux)
	Wise2_flux_err_Vals.append(rand_err_for_W2_flux)
	Wise3_flux_err_Vals.append(rand_err_for_W3_flux)
	Wise4_flux_err_Vals.append(rand_err_for_W4_flux)
	Johnson_B_flux_err_Vals.append(rand_err_for_Johnson_B_flux)
	Johnson_V_flux_err_Vals.append(rand_err_for_Johnson_V_flux)
	FUV_flux_err_Vals.append(rand_err_for_FUV_flux)
	NUV_flux_err_Vals.append(rand_err_for_NUV_flux)
	J2MASS_flux_err_Vals.append(rand_err_for_J_flux_2MASS)
	H2MASS_flux_err_Vals.append(rand_err_for_H_flux_2MASS)
	K2MASS_flux_err_Vals.append(rand_err_for_K_flux_2MASS)



temp_parameters = list(map(list, zip(Temperature_values,Gravity_values,EBV_values,Rv_values,Distance_Pc_values,sdss_u_flux_Vals,sdss_g_flux_Vals,sdss_r_flux_Vals,sdss_i_flux_Vals,sdss_z_flux_Vals,Wise1_flux_Vals,Wise2_flux_Vals,Wise3_flux_Vals,Wise4_flux_Vals,Johnson_B_flux_Vals,Johnson_V_flux_Vals,FUV_flux_Vals,NUV_flux_Vals,J2MASS_flux_Vals,H2MASS_flux_Vals,K2MASS_flux_Vals)))
parameters = tuple(temp_parameters)


#def f(TEMP,GRAV,EBV,RV,DISTANCE,J2MASS,H2MASS,K2MASS,sdssu,sdssg,sdssr,sdssi,sdssz,JohnsonB,JohnsonV,FUV,NUV,WISE1,WISE2,WISE3,WISE4):
def f(params):
	#for iteration_number in range(0,len(TEMP)):

				#SELECTING MODELS FOR 4_CLOSEST_MODEL_INTERPOLATION

		for model in numbered_model_list:				#Just need to rework the choosing of 4 closest models with modesl in /media/intern2/DD82-B92F/RESEARCH_FALL_2022/stitched_blackbody_binned_models/
			model_temp = float(model[0]+model[1]+model[2])		#numbered_model_list is 'temp1temp2temp3grav1grav2grav3'
			compare_temp = 100*model_temp
			if compare_temp >= params[0]:
				break

		for model in numbered_model_list:
			model_grav = float(model[3]+model[4]+model[5])
			compare_grav = model_grav/100
			if compare_grav >= params[1]:
				break

		max_temp = model_temp
		max_grav = model_grav

		if max_temp == 160:
			temp_low = 16000.0
			temp_high = 16000.0
			t_low = str(160)
			t_high = str(160)
			if max_grav==400:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.00
				grav_low_thigh = 4.00
				grav_high_thigh = 4.00
				glow_tlow = str(400)
				ghigh_tlow = str(400)
				glow_thigh = str(400)
				ghigh_thigh = str(400)
			if max_grav==425:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.25
				grav_low_thigh = 4.00
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(425)
				glow_thigh = str(400)
				ghigh_thigh = str(425)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.25
				grav_high_thigh = 4.50
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(425)
				ghigh_thigh = str(450)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.50
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(450)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)

		if max_temp == 180:
			temp_low = 16000.0
			temp_high = 18000.0
			t_low = str(160)
			t_high = str(180)
			if max_grav==400:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.00
				grav_low_thigh = 4.00
				grav_high_thigh = 4.00
				glow_tlow = str(400)
				ghigh_tlow = str(400)
				glow_thigh = str(400)
				ghigh_thigh = str(400)
			if max_grav==425:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.25
				grav_low_thigh = 4.00
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(425)
				glow_thigh = str(400)
				ghigh_thigh = str(425)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.25
				grav_high_thigh = 4.50
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(425)
				ghigh_thigh = str(450)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.50
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(450)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)

		if max_temp == 200:
			temp_low = 18000.0
			temp_high = 20000.0
			t_low = str(180)
			t_high = str(200)
			if max_grav==400:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.00
				grav_low_thigh = 4.00
				grav_high_thigh = 4.00
				glow_tlow = str(400)
				ghigh_tlow = str(400)
				glow_thigh = str(400)
				ghigh_thigh = str(400)
			if max_grav==425:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.25
				grav_low_thigh = 4.00
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(425)
				glow_thigh = str(400)
				ghigh_thigh = str(425)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.25
				grav_high_thigh = 4.50
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(425)
				ghigh_thigh = str(450)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.50
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(450)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.75
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(575)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.25
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(525)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)


		if max_temp == 220:
			temp_low = 20000.0
			temp_high = 22000.0
			t_low = str(200)
			t_high = str(220)
			if max_grav==400:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.00
				grav_low_thigh = 4.00
				grav_high_thigh = 4.00
				glow_tlow = str(400)
				ghigh_tlow = str(400)
				glow_thigh = str(400)
				ghigh_thigh = str(400)
			if max_grav==425:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.25
				grav_low_thigh = 4.00
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(425)
				glow_thigh = str(400)
				ghigh_thigh = str(425)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.25
				grav_high_thigh = 4.50
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(425)
				ghigh_thigh = str(450)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.50
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(450)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.75
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(575)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(525)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)

		if max_temp == 240:
			temp_low = 22000.0
			temp_high = 24000.0
			t_low = str(220)
			t_high = str(240)
			if max_grav==400:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.00
				grav_low_thigh = 4.00
				grav_high_thigh = 4.00
				glow_tlow = str(400)
				ghigh_tlow = str(400)
				glow_thigh = str(400)
				ghigh_thigh = str(400)
			if max_grav==425:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.25
				grav_low_thigh = 4.00
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(425)
				glow_thigh = str(400)
				ghigh_thigh = str(425)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.25
				grav_high_thigh = 4.50
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(425)
				ghigh_thigh = str(450)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.50
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(450)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)


		if max_temp == 260:
			temp_low = 24000.0
			temp_high = 26000.0
			t_low = str(240)
			t_high = str(260)
			if max_grav==400:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.00
				grav_low_thigh = 4.00
				grav_high_thigh = 4.00
				glow_tlow = str(400)
				ghigh_tlow = str(400)
				glow_thigh = str(400)
				ghigh_thigh = str(400)
			if max_grav==425:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.25
				grav_low_thigh = 4.00
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(425)
				glow_thigh = str(400)
				ghigh_thigh = str(425)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.25
				grav_high_thigh = 4.50
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(425)
				ghigh_thigh = str(450)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.50
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(450)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)


		if max_temp == 280:
			temp_low = 26000.0
			temp_high = 28000.0
			t_low = str(260)
			t_high = str(280)
			if max_grav==400:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.00
				grav_low_thigh = 4.00
				grav_high_thigh = 4.00
				glow_tlow = str(400)
				ghigh_tlow = str(400)
				glow_thigh = str(400)
				ghigh_thigh = str(400)
			if max_grav==425:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.25
				grav_low_thigh = 4.00
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(425)
				glow_thigh = str(400)
				ghigh_thigh = str(425)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.25
				grav_high_thigh = 4.50
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(425)
				ghigh_thigh = str(450)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.50
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(450)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)


		if max_temp == 300:
			temp_low = 28000.0
			temp_high = 30000.0
			t_low = str(280)
			t_high = str(300)
			if max_grav==400:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.00
				grav_low_thigh = 4.00
				grav_high_thigh = 4.00
				glow_tlow = str(400)
				ghigh_tlow = str(400)
				glow_thigh = str(400)
				ghigh_thigh = str(400)
			if max_grav==425:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.25
				grav_low_thigh = 4.00
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(425)
				glow_thigh = str(400)
				ghigh_thigh = str(425)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.25
				grav_high_thigh = 4.50
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(425)
				ghigh_thigh = str(450)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.50
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(450)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)

		if max_temp == 320:
			temp_low = 30000.0
			temp_high = 32000.0
			t_low = str(300)
			t_high = str(320)
			if max_grav==400:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.00
				grav_low_thigh = 4.00
				grav_high_thigh = 4.00
				glow_tlow = str(400)
				ghigh_tlow = str(400)
				glow_thigh = str(400)
				ghigh_thigh = str(400)
			if max_grav==425:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.25
				grav_low_thigh = 4.00
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(425)
				glow_thigh = str(400)
				ghigh_thigh = str(425)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.25
				grav_high_thigh = 4.50
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(425)
				ghigh_thigh = str(450)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.50
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(450)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)



		if max_temp == 340:
			temp_low = 32000.0
			temp_high = 34000.0
			t_low = str(320)
			t_high = str(340)
			if max_grav==400:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.00
				grav_low_thigh = 4.00
				grav_high_thigh = 4.00
				glow_tlow = str(400)
				ghigh_tlow = str(400)
				glow_thigh = str(400)
				ghigh_thigh = str(400)
			if max_grav==425:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.25
				grav_low_thigh = 4.00
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(425)
				glow_thigh = str(400)
				ghigh_thigh = str(425)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.25
				grav_high_thigh = 4.50
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(425)
				ghigh_thigh = str(450)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.50
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(450)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)


		if max_temp == 360:
			temp_low = 34000.0
			temp_high = 36000.0
			t_low = str(340)
			t_high = str(360)
			if max_grav==400:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.00
				grav_low_thigh = 4.00
				grav_high_thigh = 4.00
				glow_tlow = str(400)
				ghigh_tlow = str(400)
				glow_thigh = str(400)
				ghigh_thigh = str(400)
			if max_grav==425:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.25
				grav_low_thigh = 4.00
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(425)
				glow_thigh = str(400)
				ghigh_thigh = str(425)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.25
				grav_high_thigh = 4.50
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(425)
				ghigh_thigh = str(450)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.50
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(450)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)

		if max_temp == 380:
			temp_low = 36000.0
			temp_high = 38000.0
			t_low = str(360)
			t_high = str(380)
			if max_grav==400:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.00
				grav_low_thigh = 4.00
				grav_high_thigh = 4.00
				glow_tlow = str(400)
				ghigh_tlow = str(400)
				glow_thigh = str(400)
				ghigh_thigh = str(400)
			if max_grav==425:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.25
				grav_low_thigh = 4.00
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(425)
				glow_thigh = str(400)
				ghigh_thigh = str(425)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.25
				grav_high_thigh = 4.50
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(425)
				ghigh_thigh = str(450)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.50
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(450)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)

		if max_temp == 400:
			temp_low = 38000.0
			temp_high = 40000.0
			t_low = str(380)
			t_high = str(400)
			if max_grav==400:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.00
				grav_low_thigh = 4.00
				grav_high_thigh = 4.00
				glow_tlow = str(400)
				ghigh_tlow = str(400)
				glow_thigh = str(400)
				ghigh_thigh = str(400)
			if max_grav==425:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.25
				grav_low_thigh = 4.00
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(425)
				glow_thigh = str(400)
				ghigh_thigh = str(425)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.25
				grav_high_thigh = 4.50
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(425)
				ghigh_thigh = str(450)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.50
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(450)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)


		if max_temp == 450:
			temp_low = 40000.0
			temp_high = 45000.0
			t_low = str(400)
			t_high = str(450)
			if max_grav==400:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.00
				grav_low_thigh = 4.25
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(400)
				glow_thigh = str(425)
				ghigh_thigh = str(425)
			if max_grav==425:
				grav_low_tlow = 4.00
				grav_high_tlow = 4.25
				grav_low_thigh = 4.25
				grav_high_thigh = 4.25
				glow_tlow = str(400)
				ghigh_tlow = str(425)
				glow_thigh = str(425)
				ghigh_thigh = str(425)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.25
				grav_high_thigh = 4.50
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(425)
				ghigh_thigh = str(450)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.50
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(450)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)


		if max_temp == 500:
			temp_low = 45000.0
			temp_high = 50000.0
			t_low = str(450)
			t_high = str(500)
			if max_grav==400:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.25
				grav_low_thigh = 4.75
				grav_high_thigh = 4.75
				glow_tlow = str(425)
				ghigh_tlow = str(425)
				glow_thigh = str(475)
				ghigh_thigh = str(475)
			if max_grav==425:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.25
				grav_low_thigh = 4.75
				grav_high_thigh = 4.75
				glow_tlow = str(425)
				ghigh_tlow = str(425)
				glow_thigh = str(475)
				ghigh_thigh = str(475)
			if max_grav==450:
				grav_low_tlow = 4.25
				grav_high_tlow = 4.50
				grav_low_thigh = 4.75
				grav_high_thigh = 4.75
				glow_tlow = str(425)
				ghigh_tlow = str(450)
				glow_thigh = str(475)
				ghigh_thigh = str(475)
			if max_grav==475:
				grav_low_tlow = 4.50
				grav_high_tlow = 4.75
				grav_low_thigh = 4.75
				grav_high_thigh = 4.75
				glow_tlow = str(450)
				ghigh_tlow = str(475)
				glow_thigh = str(475)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)


		if Teff >= 50000.0:
			temp_low = 50000.0
			temp_high = 50000.0
			t_low = str(500)
			t_high = str(500)
			if max_grav==400:
				grav_low_tlow = 4.75
				grav_high_tlow = 4.75
				grav_low_thigh = 4.75
				grav_high_thigh = 4.75
				glow_tlow = str(475)
				ghigh_tlow = str(475)
				glow_thigh = str(475)
				ghigh_thigh = str(475)
			if max_grav==425:
				grav_low_tlow = 4.75
				grav_high_tlow = 4.75
				grav_low_thigh = 4.75
				grav_high_thigh = 4.75
				glow_tlow = str(475)
				ghigh_tlow = str(475)
				glow_thigh = str(475)
				ghigh_thigh = str(475)
			if max_grav==450:
				grav_low_tlow = 4.75
				grav_high_tlow = 4.75
				grav_low_thigh = 4.75
				grav_high_thigh = 4.75
				glow_tlow = str(475)
				ghigh_tlow = str(475)
				glow_thigh = str(475)
				ghigh_thigh = str(475)
			if max_grav==475:
				grav_low_tlow = 4.75
				grav_high_tlow = 4.75
				grav_low_thigh = 4.75
				grav_high_thigh = 4.75
				glow_tlow = str(475)
				ghigh_tlow = str(475)
				glow_thigh = str(475)
				ghigh_thigh = str(475)
			if max_grav==500:
				grav_low_tlow = 4.75
				grav_high_tlow = 5.00
				grav_low_thigh = 4.75
				grav_high_thigh = 5.00
				glow_tlow = str(475)
				ghigh_tlow = str(500)
				glow_thigh = str(475)
				ghigh_thigh = str(500)
			if max_grav==525:
				grav_low_tlow = 5.00
				grav_high_tlow = 5.25
				grav_low_thigh = 5.00
				grav_high_thigh = 5.25
				glow_tlow = str(500)
				ghigh_tlow = str(525)
				glow_thigh = str(500)
				ghigh_thigh = str(525)
			if max_grav==550:
				grav_low_tlow = 5.25
				grav_high_tlow = 5.50
				grav_low_thigh = 5.25
				grav_high_thigh = 5.50
				glow_tlow = str(525)
				ghigh_tlow = str(550)
				glow_thigh = str(525)
				ghigh_thigh = str(550)
			if max_grav==575:
				grav_low_tlow = 5.50
				grav_high_tlow = 5.75
				grav_low_thigh = 5.50
				grav_high_thigh = 5.75
				glow_tlow = str(550)
				ghigh_tlow = str(575)
				glow_thigh = str(550)
				ghigh_thigh = str(575)
			if max_grav==600:
				grav_low_tlow = 5.75
				grav_high_tlow = 6.00
				grav_low_thigh = 5.75
				grav_high_thigh = 6.00
				glow_tlow = str(575)
				ghigh_tlow = str(600)
				glow_thigh = str(575)
				ghigh_thigh = str(600)
			if max_grav==625:
				grav_low_tlow = 6.00
				grav_high_tlow = 6.25
				grav_low_thigh = 6.00
				grav_high_thigh = 6.25
				glow_tlow = str(600)
				ghigh_tlow = str(625)
				glow_thigh = str(600)
				ghigh_thigh = str(625)
			if compare_grav<=LOGg:
				grav_low_tlow = 6.25
				grav_high_tlow = 6.25
				grav_low_thigh = 6.25
				grav_high_thigh = 6.25
				glow_tlow = str(625)
				ghigh_tlow = str(625)
				glow_thigh = str(625)
				ghigh_thigh = str(625)



		models_to_convolve = []
		if temp_low == temp_high and grav_low_tlow == grav_high_tlow:
			model_one_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			models_to_convolve.append(model_one_to_convolve)
		if temp_low == temp_high and grav_low_tlow != grav_high_tlow:
			model_one_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			model_two_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			models_to_convolve.append(model_one_to_convolve)
			models_to_convolve.append(model_two_to_convolve)
		if temp_low != temp_high and grav_low_tlow == grav_high_tlow and grav_low_thigh == grav_high_thigh:
			model_one_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			model_two_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			models_to_convolve.append(model_one_to_convolve)
			models_to_convolve.append(model_two_to_convolve)
		if temp_low != temp_high and grav_low_tlow == grav_high_tlow and grav_low_thigh != grav_high_thigh:
			model_one_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			model_two_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			model_three_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			models_to_convolve.append(model_one_to_convolve)
			models_to_convolve.append(model_two_to_convolve)
			models_to_convolve.append(model_three_to_convolve)
		if temp_low != temp_high and grav_low_tlow != grav_high_tlow and grav_low_thigh == grav_high_thigh:
			model_one_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			model_two_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			model_three_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			models_to_convolve.append(model_one_to_convolve)
			models_to_convolve.append(model_two_to_convolve)
			models_to_convolve.append(model_three_to_convolve)
		if temp_low != temp_high and grav_low_tlow != grav_high_tlow and grav_low_thigh != grav_high_thigh:
			model_one_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			model_two_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			model_three_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			model_four_to_convolve = 'binned_stitched_t'+t_low+'g'+glow_tlow+'sdbv00.txt'
			models_to_convolve.append(model_one_to_convolve)
			models_to_convolve.append(model_two_to_convolve)
			models_to_convolve.append(model_three_to_convolve)
			models_to_convolve.append(model_four_to_convolve)
			#MODEL CONVOLUTION
		if len(models_to_convolve) == 1:
			convolved_model_spectrum = np.loadtxt(model_one_to_convolve)
			#print(convolved_model_spectrum)
		if len(models_to_convolve) == 2:
			a = np.loadtxt(models_to_convolve[0])
			b = np.loadtxt(models_to_convolve[1])
			a_flux = a[:,1]
			b_flux = b[:,1]
			convolved_lambda = a[:,0]
			convolved_flux = []
			if temp_low == temp_high:
				convolve_a = np.sqrt(((Teff-temp_low)**2)+((LOGg-grav_low_tlow)**2))
				convolve_b = np.sqrt(((Teff-temp_low)**2)+((LOGg-grav_high_tlow)**2))
				weight_a = (1/convolve_a)/((1/convolve_a)+(1/convolve_b))
				weight_b = (1/convolve_b)/((1/convolve_a)+(1/convolve_b))
				new_a_flux = [item*weight_a for item in a_flux]
				new_b_flux = [item*weight_b for item in b_flux]
				for lam in range(0,len(new_a_flux)):
					convolved_flux.append(new_a_flux[lam]+new_b_flux[lam])
				convolved_model_spectrum = np.array([convolved_lambda,convolved_flux]).T
				#print(convolved_model_spectrum)
			if temp_low != temp_high:
				convolve_a = np.sqrt(((Teff-temp_low)**2)+((LOGg-grav_low_tlow)**2))
				convolve_b = np.sqrt(((Teff-temp_high)**2)+((LOGg-grav_low_thigh)**2))
				weight_a = (1/convolve_a)/((1/convolve_a)+(1/convolve_b))
				weight_b = (1/convolve_b)/((1/convolve_a)+(1/convolve_b))
				new_a_flux = [item*weight_a for item in a_flux]
				new_b_flux = [item*weight_b for item in b_flux]
				for lam in range(0,len(new_a_flux)):
					convolved_flux.append(new_a_flux[lam]+new_b_flux[lam])
				convolved_model_spectrum = np.array([convolved_lambda,convolved_flux]).T
				#print(convolved_model_spectrum)
		if len(models_to_convolve) == 3:
			a = np.loadtxt(models_to_convolve[0])
			b = np.loadtxt(models_to_convolve[1])
			c = np.loadtxt(models_to_convolve[2])
			a_flux = a[:,1]
			b_flux = b[:,1]
			c_flux = c[:,1]
			convolved_lambda = a[:,0]
			convolved_flux = []
			if temp_low != temp_high and grav_low_tlow == grav_high_tlow and grav_low_thigh != grav_high_thigh:
				convolve_a = np.sqrt(((Teff-temp_low)**2)+((LOGg-grav_low_tlow)**2))
				convolve_b = np.sqrt(((Teff-temp_high)**2)+((LOGg-grav_low_thigh)**2))
				convolve_c = np.sqrt(((Teff-temp_high)**2)+((LOGg-grav_high_thigh)**2))
				weight_total = (1/convolve_a)+(1/convolve_b)+(1/convolve_c)
				weight_a = (1/convolve_a)/weight_total
				weight_b = (1/convolve_b)/weight_total
				weight_c = (1/convolve_c)/weight_total
				new_a_flux = [item*weight_a for item in a_flux]
				new_b_flux = [item*weight_b for item in b_flux]
				new_c_flux = [item*weight_c for item in c_flux]
				for lam in range(0,len(new_a_flux)):
					convolved_flux.append(new_a_flux[lam]+new_b_flux[lam]+new_c_flux[lam])
				convolved_model_spectrum = np.array([convolved_lambda,convolved_flux]).T
				#print(convolved_model_spectrum)
			if temp_low != temp_high and grav_low_tlow != grav_high_tlow and grav_low_thigh == grav_high_thigh:
				convolve_a = np.sqrt(((Teff-temp_low)**2)+((LOGg-grav_low_tlow)**2))
				convolve_b = np.sqrt(((Teff-temp_low)**2)+((LOGg-grav_high_tlow)**2))
				convolve_b = np.sqrt(((Teff-temp_high)**2)+((LOGg-grav_low_thigh)**2))
				weight_total = (1/convolve_a)+(1/convolve_b)+(1/convolve_c)
				weight_a = (1/convolve_a)/weight_total
				weight_b = (1/convolve_b)/weight_total
				weight_c = (1/convolve_c)/weight_total
				new_a_flux = [item*weight_a for item in a_flux]
				new_b_flux = [item*weight_b for item in b_flux]
				new_c_flux = [item*weight_c for item in c_flux]
				for lam in range(0,len(new_a_flux)):
					convolved_flux.append(new_a_flux[lam]+new_b_flux[lam]+new_c_flux[lam])
				convolved_model_spectrum = np.array([convolved_lambda,convolved_flux]).T
				#print(convolved_model_spectrum)
		if len(models_to_convolve) == 4:
			a = np.loadtxt(models_to_convolve[0])
			b = np.loadtxt(models_to_convolve[1])
			c = np.loadtxt(models_to_convolve[2])
			d = np.loadtxt(models_to_convolve[3])
			a_flux = a[:,1]
			b_flux = b[:,1]
			c_flux = c[:,1]
			d_flux = d[:,1]
			convolved_lambda = a[:,0]
			convolved_flux = []
			convolve_a = np.sqrt(((Teff-temp_low)**2)+((LOGg-grav_low_tlow)**2))
			convolve_b = np.sqrt(((Teff-temp_low)**2)+((LOGg-grav_high_tlow)**2))
			convolve_c = np.sqrt(((Teff-temp_high)**2)+((LOGg-grav_low_thigh)**2))
			convolve_d = np.sqrt(((Teff-temp_high)**2)+((LOGg-grav_high_thigh)**2))
			weight_total = (1/convolve_a)+(1/convolve_b)+(1/convolve_c)+(1/convolve_d)
			weight_a = (1/convolve_a)/weight_total
			weight_b = (1/convolve_b)/weight_total
			weight_c = (1/convolve_c)/weight_total
			weight_d = (1/convolve_d)/weight_total
			new_a_flux = [item*weight_a for item in a_flux]
			new_b_flux = [item*weight_b for item in b_flux]
			new_c_flux = [item*weight_c for item in c_flux]
			new_d_flux = [item*weight_d for item in d_flux]
			for lam in range(0,len(new_a_flux)):
				convolved_flux.append(new_a_flux[lam]+new_b_flux[lam]+new_c_flux[lam]+new_d_flux[lam])
			convolved_model_spectrum = np.array([convolved_lambda,convolved_flux]).T
			#print(convolved_model_spectrum)

				#Redden Convolved Model
		redx = []
		redA = []
		redfactor = []
		redcount = 0
		redlambdas = convolved_model_spectrum[:,0]
		redfluxes = convolved_model_spectrum[:,1]
		for redwl in redlambdas:
			redx.append(1/(redwl/10000))
			#if redwl == 0.0:
				#continue
				#redcount = redcount + 1
			if redx[redcount] <= 0.3:
				redA.append('out of range')
				redfactor.append(1)
				redcount = redcount + 1
			elif redx[redcount] >= 0.3 and redx[redcount] <= 1.1:
				reda,redb = 0.574*redx[redcount]**1.61,-0.527*redx[redcount]**1.61
				redA.append(params[2]*params[3]*(reda+redb/params[3]))
				redfactor.append(10**(-0.4*redA[redcount]))
				redcount = redcount + 1
			elif redx[redcount] >= 1.1 and redx[redcount] <= 3.3:
				redy = redx[redcount]-1.82
				reda,redb = 1.0+0.17699*redy-0.50447*redy**2-0.02427*redy**3+0.72085*redy**4+0.050447*redy**5-0.77530*redy**6+0.32999*redy**7,1.41338*redy+2.28305*redy**2+1.07266*redy**3-5.38434*redy**4-0.62251*redy**5+5.30260*redy**6-2.09002*redy**7
				redA.append(params[2]*params[3]*(reda+redb/params[3]))
				redfactor.append(10**(-0.4*redA[redcount]))
				redcount = redcount + 1
			elif redx[redcount] >= 3.3 and redx[redcount] <= 8.0:
				if redx[redcount] < 5.9:
					redFa,redFb = 0.0,0.0
				else:
					redFa,redFb = -0.04473*(redx[redcount]-5.9)**2-0.009779*(redx[redcount]-5.9)**3,0.21300*(redx[redcount]-5.9)**2+0.120700*(redx[redcount]-5.9)**3
				reda,redb = 1.752-0.316*redx[redcount]-0.104/((redx[redcount]-4.67)**2+0.341)+redFa,-3.090+1.825*redx[redcount]+1.206/((redx[redcount]-4.62)**2+0.263)+redFb
				redA.append(params[2]*params[3]*(reda+redb/params[3]))
				redfactor.append(10**(-0.4*redA[redcount]))
				redcount = redcount + 1
			elif redx[redcount] >= 8.0 and redx[redcount] <= 10.0:
				reda,redb = -1.073-0.628*(redx[redcount]-8.0)+0.137*(redx[redcount]-8.0)**2-0.070*(redx[redcount]-8.0)**3,13.670+4.257*(redx[redcount]-8.0)-0.420*(redx[redcount]-8.0)**2+0.374*(redx[redcount]-8.0)**3
				redA.append(params[2]*params[3]*(reda+redb/params[3]))
				redfactor.append(10**(-0.4*redA[redcount]))
				redcount = redcount + 1
			elif redx[redcount] >= 10.0:
				redA.append('out of range')
				redfactor.append(1)
				redcount = redcount + 1


		reddened_flux,redcount = [],0
		for redf in redfluxes:
			reddened_flux.append(redf*redfactor[redcount])
			redcount = redcount + 1
		reddened_convolved_model_spectrum = array([redlambdas,reddened_flux]).T


			#SMOOTH REDDENED_CONVOLVED_MODEL_SPECTRUM
		#fwhm = 12 if you go based on the narrowest absorbtion line in GD934 spectrum... Dr. Reed said to try a constant 5 Angstrom seperation
		#sigma = fwhm/sqrt(8*log(2))
		sigma = 2.1			#this corresponds to a fwhm of 5 angstroms

		smoothed_reddened_convolved_model_fluxes = gaussian_filter(reddened_convolved_model_spectrum[:,1], sigma)
		smoothed_reddened_convolved_model_spectrum = np.array([reddened_convolved_model_spectrum[:,0],smoothed_reddened_convolved_model_fluxes]).T



			#SYNTHETIC PHOTOMETRY
		smoothed_reddened_convolved_model_lams = smoothed_reddened_convolved_model_spectrum[:,0] * unit['AA']
		smoothed_reddened_convolved_model_fluxes = smoothed_reddened_convolved_model_spectrum[:,1] * unit['erg/s/cm**2/AA']
		interpolation_of_smoothed_reddened_convolved_model = interp1d(smoothed_reddened_convolved_model_spectrum[:,0],smoothed_reddened_convolved_model_spectrum[:,1],'linear')
		flux_for_2MASS_K = []
		lam_for_2MASS_K = []
		temp_lams_for_2MASS_K = np.linspace(1.939000000000000000e+04,2.384000000000000000e+04,4451)
		for lam in temp_lams_for_2MASS_K:
			flux_for_2MASS_K.append(interpolation_of_smoothed_reddened_convolved_model(lam))
			lam_for_2MASS_K.append(lam)
		array_for_2MASS_K = np.array([lam_for_2MASS_K,flux_for_2MASS_K]).T
		lams_for_2MASS_K = array_for_2MASS_K[:,0] * unit['AA']
		fluxes_for_2MASS_K = array_for_2MASS_K[:,1] * unit['erg/s/cm**2/AA']

		#gaia_g_quantity = gaia_g_filter.getFlux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
		#gaia_g_temp_str = f"{gaia_g_quantity}"
		#gaia_g_string = gaia_g_temp_str[0]+gaia_g_temp_str[1]+gaia_g_temp_str[2]+gaia_g_temp_str[3]+gaia_g_temp_str[4]+gaia_g_temp_str[5]+gaia_g_temp_str[6]+gaia_g_temp_str[7]+gaia_g_temp_str[8]+gaia_g_temp_str[9]+gaia_g_temp_str[10]+gaia_g_temp_str[11]+gaia_g_temp_str[12]+gaia_g_temp_str[13]
		#gaia_g_synthetic_flux_erg = np.float64(gaia_g_string)

		if type(params[18]) != str:
			J_2MASS_quantity = J_2MASS_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			J_2MASS_temp_str = f"{J_2MASS_quantity}"
			J_2MASS_string = J_2MASS_temp_str[0]+J_2MASS_temp_str[1]+J_2MASS_temp_str[2]+J_2MASS_temp_str[3]+J_2MASS_temp_str[4]+J_2MASS_temp_str[5]+J_2MASS_temp_str[6]+J_2MASS_temp_str[7]+J_2MASS_temp_str[8]+J_2MASS_temp_str[9]+J_2MASS_temp_str[10]+J_2MASS_temp_str[11]+J_2MASS_temp_str[12]
			J_2MASS_real_string = []
			for character in range(0,len(J_2MASS_string)):
				if J_2MASS_string[character] != ' ' and J_2MASS_string[character] != '"':
					J_2MASS_real_string.append(J_2MASS_string[character])
				else:
					break
			J_2MASS_realL_string = J_2MASS_real_string[0]+J_2MASS_real_string[1]+J_2MASS_real_string[2]+J_2MASS_real_string[3]+J_2MASS_real_string[4]+J_2MASS_real_string[5]+J_2MASS_real_string[6]+J_2MASS_real_string[7]+J_2MASS_real_string[8]+J_2MASS_real_string[9]+J_2MASS_real_string[10]+J_2MASS_real_string[11]+J_2MASS_real_string[12]
			J_2MASS_synthetic_flux_erg = np.float64(J_2MASS_realL_string)

		if type(params[18]) == str:
			J_2MASS_synthetic_flux_erg = 'Not Needed'

		if type(params[19]) != str:
			H_2MASS_quantity = H_2MASS_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			H_2MASS_temp_str = f"{H_2MASS_quantity}"
			H_2MASS_string = H_2MASS_temp_str[0]+H_2MASS_temp_str[1]+H_2MASS_temp_str[2]+H_2MASS_temp_str[3]+H_2MASS_temp_str[4]+H_2MASS_temp_str[5]+H_2MASS_temp_str[6]+H_2MASS_temp_str[7]+H_2MASS_temp_str[8]+H_2MASS_temp_str[9]+H_2MASS_temp_str[10]+H_2MASS_temp_str[11]+H_2MASS_temp_str[12]
			H_2MASS_real_string = []
			for character in range(0,len(H_2MASS_string)):
				if H_2MASS_string[character] != ' ' and H_2MASS_string[character] != '"':
					H_2MASS_real_string.append(H_2MASS_string[character])
				else:
					break
			H_2MASS_realL_string = H_2MASS_real_string[0]+H_2MASS_real_string[1]+H_2MASS_real_string[2]+H_2MASS_real_string[3]+H_2MASS_real_string[4]+H_2MASS_real_string[5]+H_2MASS_real_string[6]+H_2MASS_real_string[7]+H_2MASS_real_string[8]+H_2MASS_real_string[9]+H_2MASS_real_string[10]+H_2MASS_real_string[11]+H_2MASS_real_string[12]
			H_2MASS_synthetic_flux_erg = np.float64(H_2MASS_realL_string)

		if type(params[19]) == str:
			H_2MASS_synthetic_flux_erg = 'Not Needed'

		if type(params[20]) != str:
			K_2MASS_quantity = K_2MASS_filter.get_flux(lams_for_2MASS_K,fluxes_for_2MASS_K)
			K_2MASS_temp_str = f"{K_2MASS_quantity}"
			K_2MASS_string = K_2MASS_temp_str[0]+K_2MASS_temp_str[1]+K_2MASS_temp_str[2]+K_2MASS_temp_str[3]+K_2MASS_temp_str[4]+K_2MASS_temp_str[5]+K_2MASS_temp_str[6]+K_2MASS_temp_str[7]+K_2MASS_temp_str[8]+K_2MASS_temp_str[9]+K_2MASS_temp_str[10]+K_2MASS_temp_str[11]+K_2MASS_temp_str[12]
			K_2MASS_real_string = []
			for character in range(0,len(K_2MASS_string)):
				if K_2MASS_string[character] != ' ' and K_2MASS_string[character] != '"':
					K_2MASS_real_string.append(K_2MASS_string[character])
				else:
					break
			K_2MASS_realL_string = K_2MASS_real_string[0]+K_2MASS_real_string[1]+K_2MASS_real_string[2]+K_2MASS_real_string[3]+K_2MASS_real_string[4]+K_2MASS_real_string[5]+K_2MASS_real_string[6]+K_2MASS_real_string[7]+K_2MASS_real_string[8]+K_2MASS_real_string[9]+K_2MASS_real_string[10]+K_2MASS_real_string[11]+K_2MASS_real_string[12]
			K_2MASS_synthetic_flux_erg = np.float64(K_2MASS_realL_string)

		if type(params[20]) == str:
			K_2MASS_synthetic_flux_erg = 'Not Needed'

		#Johnson_U_quantity = Johnson_U_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
		#Johnson_U_temp_str = f"{Johnson_U_quantity}"
		#Johnson_U_string = Johnson_U_temp_str[0]+Johnson_U_temp_str[1]+Johnson_U_temp_str[2]+Johnson_U_temp_str[3]+Johnson_U_temp_str[4]+Johnson_U_temp_str[5]+Johnson_U_temp_str[6]+Johnson_U_temp_str[7]+Johnson_U_temp_str[8]+Johnson_U_temp_str[9]+Johnson_U_temp_str[10]+Johnson_U_temp_str[11]+Johnson_U_temp_str[12]+Johnson_U_temp_str[13]
		#Johsnon_U_synthetic_flux_erg = np.float64(Johnson_U_string)

		if type(params[14]) != str:
			Johnson_B_quantity = Johnson_B_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			Johnson_B_temp_str = f"{Johnson_B_quantity}"
			Johnson_B_string = Johnson_B_temp_str[0]+Johnson_B_temp_str[1]+Johnson_B_temp_str[2]+Johnson_B_temp_str[3]+Johnson_B_temp_str[4]+Johnson_B_temp_str[5]+Johnson_B_temp_str[6]+Johnson_B_temp_str[7]+Johnson_B_temp_str[8]+Johnson_B_temp_str[9]+Johnson_B_temp_str[10]+Johnson_B_temp_str[11]+Johnson_B_temp_str[12]
			Johnson_B_real_string = []
			for character in range (0,len(Johnson_B_string)):
				if Johnson_B_string[character] != ' ' and Johnson_B_string[character] != '"':
					Johnson_B_real_string.append(Johnson_B_string[character])
				else:
					break
			Johnson_B_realL_string = Johnson_B_real_string[0]+Johnson_B_real_string[1]+Johnson_B_real_string[2]+Johnson_B_real_string[3]+Johnson_B_real_string[4]+Johnson_B_real_string[5]+Johnson_B_real_string[6]+Johnson_B_real_string[7]+Johnson_B_real_string[8]+Johnson_B_real_string[9]+Johnson_B_real_string[10]+Johnson_B_real_string[11]+Johnson_B_real_string[12]
			Johnson_B_synthetic_flux_erg = np.float64(Johnson_B_realL_string)

		if type(params[14]) == str:
			Johnson_B_synthetic_flux_erg = 'Not Needed'

		if type(params[15]) != str:
			Johnson_V_quantity = Johnson_V_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			Johnson_V_temp_str = f"{Johnson_V_quantity}"
			Johnson_V_string = Johnson_V_temp_str[0]+Johnson_V_temp_str[1]+Johnson_V_temp_str[2]+Johnson_V_temp_str[3]+Johnson_V_temp_str[4]+Johnson_V_temp_str[5]+Johnson_V_temp_str[6]+Johnson_V_temp_str[7]+Johnson_V_temp_str[8]+Johnson_V_temp_str[9]+Johnson_V_temp_str[10]+Johnson_V_temp_str[11]+Johnson_V_temp_str[12]
			Johnson_V_real_string = []
			for character in range (0,len(Johnson_V_string)):
				if Johnson_V_string[character] != ' ' and Johnson_V_string[character] != '"':
					Johnson_V_real_string.append(Johnson_V_string[character])
				else:
					break
			Johnson_V_realL_string = Johnson_V_real_string[0]+Johnson_V_real_string[1]+Johnson_V_real_string[2]+Johnson_V_real_string[3]+Johnson_V_real_string[4]+Johnson_V_real_string[5]+Johnson_V_real_string[6]+Johnson_V_real_string[7]+Johnson_V_real_string[8]+Johnson_V_real_string[9]+Johnson_V_real_string[10]+Johnson_V_real_string[11]+Johnson_V_real_string[12]
			Johnson_V_synthetic_flux_erg = np.float64(Johnson_V_realL_string)

		if type(params[15]) == str:
			Johnson_V_synthetic_flux_erg = 'Not Needed'

		if type(params[5]) != str:
			sdss_u_quantity = sdss_u_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			sdss_u_temp_str = f"{sdss_u_quantity}"
			sdss_u_string = sdss_u_temp_str[0]+sdss_u_temp_str[1]+sdss_u_temp_str[2]+sdss_u_temp_str[3]+sdss_u_temp_str[4]+sdss_u_temp_str[5]+sdss_u_temp_str[6]+sdss_u_temp_str[7]+sdss_u_temp_str[8]+sdss_u_temp_str[9]+sdss_u_temp_str[10]+sdss_u_temp_str[11]+sdss_u_temp_str[12]
			sdss_u_real_string = []
			for character in range (0,len(sdss_u_string)):
				if sdss_u_string[character] != ' ' and sdss_u_string[character] != '"':
					sdss_u_real_string.append(sdss_u_string[character])
				else:
					break
			sdss_u_realL_string = sdss_u_real_string[0]+sdss_u_real_string[1]+sdss_u_real_string[2]+sdss_u_real_string[3]+sdss_u_real_string[4]+sdss_u_real_string[5]+sdss_u_real_string[6]+sdss_u_real_string[7]+sdss_u_real_string[8]+sdss_u_real_string[9]+sdss_u_real_string[10]+sdss_u_real_string[11]+sdss_u_real_string[12]
			sdss_u_synthetic_flux_erg = np.float64(sdss_u_realL_string)

		if type(params[5]) == str:
			sdss_u_synthetic_flux_erg = 'Not Needed'

		if type(params[6]) != str:
			sdss_g_quantity = sdss_g_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			sdss_g_temp_str = f"{sdss_g_quantity}"
			sdss_g_string = sdss_g_temp_str[0]+sdss_g_temp_str[1]+sdss_g_temp_str[2]+sdss_g_temp_str[3]+sdss_g_temp_str[4]+sdss_g_temp_str[5]+sdss_g_temp_str[6]+sdss_g_temp_str[7]+sdss_g_temp_str[8]+sdss_g_temp_str[9]+sdss_g_temp_str[10]+sdss_g_temp_str[11]+sdss_g_temp_str[12]
			sdss_g_real_string = []
			for character in range (0,len(sdss_g_string)):
				if sdss_g_string[character] != ' ' and sdss_g_string[character] != '"':
					sdss_g_real_string.append(sdss_g_string[character])
				else:
					break
			sdss_g_realL_string = sdss_g_real_string[0]+sdss_g_real_string[1]+sdss_g_real_string[2]+sdss_g_real_string[3]+sdss_g_real_string[4]+sdss_g_real_string[5]+sdss_g_real_string[6]+sdss_g_real_string[7]+sdss_g_real_string[8]+sdss_g_real_string[9]+sdss_g_real_string[10]+sdss_g_real_string[11]+sdss_g_real_string[12]
			sdss_g_synthetic_flux_erg = np.float64(sdss_g_realL_string)

		if type(params[6]) == str:
			sdss_g_synthetic_flux_erg = 'Not Needed'

		if type(params[7]) != str:
			sdss_r_quantity = sdss_r_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			sdss_r_temp_str = f"{sdss_r_quantity}"
			sdss_r_string = sdss_r_temp_str[0]+sdss_r_temp_str[1]+sdss_r_temp_str[2]+sdss_r_temp_str[3]+sdss_r_temp_str[4]+sdss_r_temp_str[5]+sdss_r_temp_str[6]+sdss_r_temp_str[7]+sdss_r_temp_str[8]+sdss_r_temp_str[9]+sdss_r_temp_str[10]+sdss_r_temp_str[11]+sdss_r_temp_str[12]
			sdss_r_real_string = []
			for character in range (0,len(sdss_r_string)):
				if sdss_r_string[character] != ' ' and sdss_r_string[character] != '"':
					sdss_r_real_string.append(sdss_r_string[character])
				else:
					break
			sdss_r_realL_string = sdss_r_real_string[0]+sdss_r_real_string[1]+sdss_r_real_string[2]+sdss_r_real_string[3]+sdss_r_real_string[4]+sdss_r_real_string[5]+sdss_r_real_string[6]+sdss_r_real_string[7]+sdss_r_real_string[8]+sdss_r_real_string[9]+sdss_r_real_string[10]+sdss_r_real_string[11]+sdss_r_real_string[12]
			sdss_r_synthetic_flux_erg = np.float64(sdss_r_realL_string)

		if type(params[7]) == str:
			sdss_r_synthetic_flux_erg = 'Not Needed'

		if type(params[8]) != str:
			sdss_i_quantity = sdss_i_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			sdss_i_temp_str = f"{sdss_i_quantity}"
			sdss_i_string = sdss_i_temp_str[0]+sdss_i_temp_str[1]+sdss_i_temp_str[2]+sdss_i_temp_str[3]+sdss_i_temp_str[4]+sdss_i_temp_str[5]+sdss_i_temp_str[6]+sdss_i_temp_str[7]+sdss_i_temp_str[8]+sdss_i_temp_str[9]+sdss_i_temp_str[10]+sdss_i_temp_str[11]+sdss_i_temp_str[12]
			sdss_i_real_string = []
			for character in range (0,len(sdss_i_string)):
				if sdss_i_string[character] != ' ' and sdss_i_string[character] != '"':
					sdss_i_real_string.append(sdss_i_string[character])
				else:
					break
			sdss_i_realL_string = sdss_i_real_string[0]+sdss_i_real_string[1]+sdss_i_real_string[2]+sdss_i_real_string[3]+sdss_i_real_string[4]+sdss_i_real_string[5]+sdss_i_real_string[6]+sdss_i_real_string[7]+sdss_i_real_string[8]+sdss_i_real_string[9]+sdss_i_real_string[10]+sdss_i_real_string[11]+sdss_i_real_string[12]
			sdss_i_synthetic_flux_erg = np.float64(sdss_i_realL_string)

		if type(params[8]) == str:
			sdss_i_synthetic_flux_erg = 'Not Needed'

		if type(params[9]) != str:
			sdss_z_quantity = sdss_z_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			sdss_z_temp_str = f"{sdss_z_quantity}"
			sdss_z_string = sdss_z_temp_str[0]+sdss_z_temp_str[1]+sdss_z_temp_str[2]+sdss_z_temp_str[3]+sdss_z_temp_str[4]+sdss_z_temp_str[5]+sdss_z_temp_str[6]+sdss_z_temp_str[7]+sdss_z_temp_str[8]+sdss_z_temp_str[9]+sdss_z_temp_str[10]+sdss_z_temp_str[11]+sdss_z_temp_str[12]
			sdss_z_real_string = []
			for character in range (0,len(sdss_z_string)):
				if sdss_z_string[character] != ' ' and sdss_z_string[character] != '"':
					sdss_z_real_string.append(sdss_z_string[character])
				else:
					break
			sdss_z_realL_string = sdss_z_real_string[0]+sdss_z_real_string[1]+sdss_z_real_string[2]+sdss_z_real_string[3]+sdss_z_real_string[4]+sdss_z_real_string[5]+sdss_z_real_string[6]+sdss_z_real_string[7]+sdss_z_real_string[8]+sdss_z_real_string[9]+sdss_z_real_string[10]+sdss_z_real_string[11]+sdss_z_real_string[12]
			sdss_z_synthetic_flux_erg = np.float64(sdss_z_realL_string)

		if type(params[9]) == str:
			sdss_z_synthetic_flux_erg = 'Not Needed'

		if type(params[10]) != str:
			WISE1_quantity = WISE1_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			WISE1_temp_str = f"{WISE1_quantity}"
			WISE1_string = WISE1_temp_str[0]+WISE1_temp_str[1]+WISE1_temp_str[2]+WISE1_temp_str[3]+WISE1_temp_str[4]+WISE1_temp_str[5]+WISE1_temp_str[6]+WISE1_temp_str[7]+WISE1_temp_str[8]+WISE1_temp_str[9]+WISE1_temp_str[10]+WISE1_temp_str[11]+WISE1_temp_str[12]
			WISE1_real_string = []
			for character in range (0,len(WISE1_string)):
				if WISE1_string[character] != ' ' and WISE1_string[character] != '"':
					WISE1_real_string.append(WISE1_string[character])
				else:
					break
			WISE1_realL_string = WISE1_real_string[0]+WISE1_real_string[1]+WISE1_real_string[2]+WISE1_real_string[3]+WISE1_real_string[4]+WISE1_real_string[5]+WISE1_real_string[6]+WISE1_real_string[7]+WISE1_real_string[8]+WISE1_real_string[9]+WISE1_real_string[10]+WISE1_real_string[11]+WISE1_real_string[12]
			WISE1_synthetic_flux_erg = np.float64(WISE1_realL_string)

		if type(params[10]) == str:
			WISE1_synthetic_flux_erg = 'Not Needed'

		if type(params[11]) != str:
			WISE2_quantity = WISE2_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			WISE2_temp_str = f"{WISE2_quantity}"
			WISE2_string = WISE2_temp_str[0]+WISE2_temp_str[1]+WISE2_temp_str[2]+WISE2_temp_str[3]+WISE2_temp_str[4]+WISE2_temp_str[5]+WISE2_temp_str[6]+WISE2_temp_str[7]+WISE2_temp_str[8]+WISE2_temp_str[9]+WISE2_temp_str[10]+WISE2_temp_str[11]+WISE2_temp_str[12]
			WISE2_real_string = []
			for character in range (0,len(WISE2_string)):
				if WISE2_string[character] != ' ' and WISE2_string[character] != '"':
					WISE2_real_string.append(WISE2_string[character])
				else:
					break
			WISE2_realL_string = WISE2_real_string[0]+WISE2_real_string[1]+WISE2_real_string[2]+WISE2_real_string[3]+WISE2_real_string[4]+WISE2_real_string[5]+WISE2_real_string[6]+WISE2_real_string[7]+WISE2_real_string[8]+WISE2_real_string[9]+WISE2_real_string[10]+WISE2_real_string[11]+WISE2_real_string[12]
			WISE2_synthetic_flux_erg = np.float64(WISE2_realL_string)

		if type(params[11]) == str:
			WISE2_synthetic_flux_erg = 'Not Needed'

		if type(params[12]) != str:
			WISE3_quantity = WISE3_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			WISE3_temp_str = f"{WISE3_quantity}"
			WISE3_string = WISE3_temp_str[0]+WISE3_temp_str[1]+WISE3_temp_str[2]+WISE3_temp_str[3]+WISE3_temp_str[4]+WISE3_temp_str[5]+WISE3_temp_str[6]+WISE3_temp_str[7]+WISE3_temp_str[8]+WISE3_temp_str[9]+WISE3_temp_str[10]+WISE3_temp_str[11]+WISE3_temp_str[12]
			WISE3_real_string = []
			for character in range (0,len(WISE3_string)):
				if WISE3_string[character] != ' ' and WISE3_string[character] != '"':
					WISE3_real_string.append(WISE3_string[character])
				else:
					break
			WISE3_realL_string = WISE3_real_string[0]+WISE3_real_string[1]+WISE3_real_string[2]+WISE3_real_string[3]+WISE3_real_string[4]+WISE3_real_string[5]+WISE3_real_string[6]+WISE3_real_string[7]+WISE3_real_string[8]+WISE3_real_string[9]+WISE3_real_string[10]+WISE3_real_string[11]+WISE3_real_string[12]
			WISE3_synthetic_flux_erg = np.float64(WISE3_realL_string)

		if type(params[12]) == str:
			WISE3_synthetic_flux_erg = 'Not Needed'

		if type(params[13]) != str:
			WISE4_quantity = WISE4_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			WISE4_temp_str = f"{WISE4_quantity}"
			WISE4_string = WISE4_temp_str[0]+WISE4_temp_str[1]+WISE4_temp_str[2]+WISE4_temp_str[3]+WISE4_temp_str[4]+WISE4_temp_str[5]+WISE4_temp_str[6]+WISE4_temp_str[7]+WISE4_temp_str[8]+WISE4_temp_str[9]+WISE4_temp_str[10]+WISE4_temp_str[11]+WISE4_temp_str[12]
			WISE4_real_string = []
			for character in range (0,len(WISE4_string)):
				if WISE4_string[character] != ' ' and WISE4_string[character] != '"':
					WISE4_real_string.append(WISE4_string[character])
				else:
					break
			WISE4_realL_string = WISE4_real_string[0]+WISE4_real_string[1]+WISE4_real_string[2]+WISE4_real_string[3]+WISE4_real_string[4]+WISE4_real_string[5]+WISE4_real_string[6]+WISE4_real_string[7]+WISE4_real_string[8]+WISE4_real_string[9]+WISE4_real_string[10]+WISE4_real_string[11]+WISE4_real_string[12]
			WISE4_synthetic_flux_erg = np.float64(WISE4_realL_string)

		if type(params[13]) == str:
			WISE4_synthetic_flux_erg = 'Not Needed'

		if type(params[16]) != str:
			fuv_galex_quantity = fuv_galex_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			fuv_galex_temp_str = f"{fuv_galex_quantity}"
			fuv_galex_string = fuv_galex_temp_str[0]+fuv_galex_temp_str[1]+fuv_galex_temp_str[2]+fuv_galex_temp_str[3]+fuv_galex_temp_str[4]+fuv_galex_temp_str[5]+fuv_galex_temp_str[6]+fuv_galex_temp_str[7]+fuv_galex_temp_str[8]+fuv_galex_temp_str[9]+fuv_galex_temp_str[10]+fuv_galex_temp_str[11]+fuv_galex_temp_str[12]
			fuv_galex_real_string = []
			for character in range (0,len(fuv_galex_string)):
				if fuv_galex_string[character] != ' ' and fuv_galex_string[character] != '"':
					fuv_galex_real_string.append(fuv_galex_string[character])
				else:
					break
			fuv_galex_realL_string = fuv_galex_real_string[0]+fuv_galex_real_string[1]+fuv_galex_real_string[2]+fuv_galex_real_string[3]+fuv_galex_real_string[4]+fuv_galex_real_string[5]+fuv_galex_real_string[6]+fuv_galex_real_string[7]+fuv_galex_real_string[8]+fuv_galex_real_string[9]+fuv_galex_real_string[10]+fuv_galex_real_string[11]+fuv_galex_real_string[12]
			fuv_galex_synthetic_flux_erg = np.float64(fuv_galex_realL_string)

		if type(params[16]) == str:
			fuv_galex_synthetic_flux_erg = 'Not Needed'

		if type(params[17]) != str:
			nuv_galex_quantity = nuv_galex_filter.get_flux(smoothed_reddened_convolved_model_lams,smoothed_reddened_convolved_model_fluxes)
			nuv_galex_temp_str = f"{nuv_galex_quantity}"
			nuv_galex_string = nuv_galex_temp_str[0]+nuv_galex_temp_str[1]+nuv_galex_temp_str[2]+nuv_galex_temp_str[3]+nuv_galex_temp_str[4]+nuv_galex_temp_str[5]+nuv_galex_temp_str[6]+nuv_galex_temp_str[7]+nuv_galex_temp_str[8]+nuv_galex_temp_str[9]+nuv_galex_temp_str[10]+nuv_galex_temp_str[11]+nuv_galex_temp_str[12]
			nuv_galex_real_string = []
			for character in range (0,len(nuv_galex_string)):
				if nuv_galex_string[character] != ' ' and nuv_galex_string[character] != '"':
					nuv_galex_real_string.append(nuv_galex_string[character])
				else:
					break
			nuv_galex_realL_string = nuv_galex_real_string[0]+nuv_galex_real_string[1]+nuv_galex_real_string[2]+nuv_galex_real_string[3]+nuv_galex_real_string[4]+nuv_galex_real_string[5]+nuv_galex_real_string[6]+nuv_galex_real_string[7]+nuv_galex_real_string[8]+nuv_galex_real_string[9]+nuv_galex_real_string[10]+nuv_galex_real_string[11]+nuv_galex_real_string[12]
			nuv_galex_synthetic_flux_erg = np.float64(nuv_galex_realL_string)

		if type(params[17]) == str:
			nuv_galex_synthetic_flux_erg = 'Not Needed'



			#(observed) Relative Errors for filters
		total_relative_err_to_sum = []
		if type(J_2MASS_synthetic_flux_erg) != str:
			rel_err_J_2MASS = vizier_J_flux_err_erg/vizier_J_flux_erg
			total_relative_err_to_sum.append(rel_err_J_2MASS)

		if type(H_2MASS_synthetic_flux_erg) != str:
			rel_err_H_2MASS = vizier_H_flux_err_erg/vizier_H_flux_erg
			total_relative_err_to_sum.append(rel_err_H_2MASS)

		if type(K_2MASS_synthetic_flux_erg) != str:
			rel_err_K_2MASS = vizier_K_flux_err_erg/vizier_K_flux_erg
			total_relative_err_to_sum.append(rel_err_H_2MASS)

		if type(Johnson_B_synthetic_flux_erg) != str:
			rel_err_Johnson_B = Johnson_B_flux_err_erg/Johnson_B_flux_erg
			total_relative_err_to_sum.append(rel_err_Johnson_B)

		if type(Johnson_V_synthetic_flux_erg) != str:
			rel_err_Johnson_V = Johnson_V_flux_err_erg/Johnson_V_flux_erg
			total_relative_err_to_sum.append(rel_err_Johnson_V)

		if type(sdss_u_synthetic_flux_erg) != str:
			rel_err_sdss_u = err_u_flux_avg_erg/u_flux_avg_erg
			total_relative_err_to_sum.append(rel_err_sdss_u)

		if type(sdss_g_synthetic_flux_erg) != str:
			rel_err_sdss_g = err_g_flux_avg_erg/g_flux_avg_erg
			total_relative_err_to_sum.append(rel_err_sdss_g)

		if type(sdss_r_synthetic_flux_erg) != str:
			rel_err_sdss_r = err_r_flux_avg_erg/r_flux_avg_erg
			total_relative_err_to_sum.append(rel_err_sdss_r)

		if type(sdss_i_synthetic_flux_erg) != str:
			rel_err_sdss_i = err_i_flux_avg_erg/i_flux_avg_erg
			total_relative_err_to_sum.append(rel_err_sdss_i)

		if type(sdss_z_synthetic_flux_erg) != str:
			rel_err_sdss_z = err_z_flux_avg_erg/z_flux_avg_erg
			total_relative_err_to_sum.append(rel_err_sdss_z)

		if type(WISE1_synthetic_flux_erg) != str:
			rel_err_WISE1 = erg_error_list[0]/erg_list[0]
			total_relative_err_to_sum.append(rel_err_WISE1)

		if type(WISE2_synthetic_flux_erg) != str:
			rel_err_WISE2 = erg_error_list[1]/erg_list[1]
			total_relative_err_to_sum.append(rel_err_WISE2)

		if type(WISE3_synthetic_flux_erg) != str:
			rel_err_WISE3 = erg_error_list[2]/erg_list[2]
			total_relative_err_to_sum.append(rel_err_WISE3)

		if type(WISE4_synthetic_flux_erg) != str:
			rel_err_WISE4 = erg_error_list[3]/erg_list[3]
			total_relative_err_to_sum.append(rel_err_WISE4)

		if type(fuv_galex_synthetic_flux_erg) != str:
			rel_err_fuv_galex = avg_galex_fuv_flux_err/avg_galex_fuv_flux
			total_relative_err_to_sum.append(rel_err_fuv_galex)

		if type(nuv_galex_synthetic_flux_erg) != str:
			rel_err_nuv_galex = avg_galex_nuv_flux_err/avg_galex_nuv_flux
			total_relative_err_to_sum.append(rel_err_nuv_galex)

		total_relative_err = np.sum(total_relative_err_to_sum)


			#FLUX RATIO / CHI SQUARED
		to_average_fr = []
		total_relative_err_to_sum = []
		if type(J_2MASS_synthetic_flux_erg) != str:
			weight_J_2MASS = 1-(rel_err_J_2MASS/total_relative_err)
			to_average_fr.append(weight_J_2MASS*(params[18]/J_2MASS_synthetic_flux_erg))

		if type(H_2MASS_synthetic_flux_erg) != str:
			weight_H_2MASS = 1-(rel_err_H_2MASS/total_relative_err)
			to_average_fr.append(weight_H_2MASS*(params[19]/H_2MASS_synthetic_flux_erg))

		if type(K_2MASS_synthetic_flux_erg) != str:
			weight_K_2MASS = 1-(rel_err_K_2MASS/total_relative_err)
			to_average_fr.append(weight_K_2MASS*(params[20]/K_2MASS_synthetic_flux_erg))

		if type(Johnson_B_synthetic_flux_erg) != str:
			weight_Johnson_B = 1-(rel_err_Johnson_B/total_relative_err)
			to_average_fr.append(weight_Johnson_B*(params[14]/Johnson_B_synthetic_flux_erg))

		if type(Johnson_V_synthetic_flux_erg) != str:
			weight_Johnson_V = 1-(rel_err_Johnson_V/total_relative_err)
			to_average_fr.append(weight_Johnson_V*(params[15]/Johnson_V_synthetic_flux_erg))

		if type(sdss_u_synthetic_flux_erg) != str:
			weight_sdss_u = 1-(rel_err_sdss_u/total_relative_err)
			to_average_fr.append(weight_sdss_u*(params[5]/sdss_u_synthetic_flux_erg))

		if type(sdss_g_synthetic_flux_erg) != str:
			weight_sdss_g = 1-(rel_err_sdss_g/total_relative_err)
			to_average_fr.append(weight_sdss_g*(params[6]/sdss_g_synthetic_flux_erg))

		if type(sdss_r_synthetic_flux_erg) != str:
			weight_sdss_r = 1-(rel_err_sdss_r/total_relative_err)
			to_average_fr.append(weight_sdss_r*(params[7]/sdss_r_synthetic_flux_erg))

		if type(sdss_i_synthetic_flux_erg) != str:
			weight_sdss_i = 1-(rel_err_sdss_i/total_relative_err)
			to_average_fr.append(weight_sdss_i*(params[8]/sdss_i_synthetic_flux_erg))

		if type(sdss_z_synthetic_flux_erg) != str:
			weight_sdss_z = 1-(rel_err_sdss_z/total_relative_err)
			to_average_fr.append(weight_sdss_z*(params[9]/sdss_z_synthetic_flux_erg))

		if type(WISE1_synthetic_flux_erg) != str:
			weight_WISE1 = 1-(rel_err_WISE1/total_relative_err)
			to_average_fr.append(weight_WISE1*(params[10]/WISE1_synthetic_flux_erg))

		if type(WISE2_synthetic_flux_erg) != str:
			weight_WISE2 = 1-(rel_err_WISE2/total_relative_err)
			to_average_fr.append(weight_WISE2*(params[11]/WISE2_synthetic_flux_erg))

		if type(WISE3_synthetic_flux_erg) != str:
			weight_WISE3 = 1-(rel_err_WISE3/total_relative_err)
			to_average_fr.append(weight_WISE3*(params[12]/WISE3_synthetic_flux_erg))

		if type(WISE4_synthetic_flux_erg) != str:
			weight_WISE4 = 1-(rel_err_WISE4/total_relative_err)
			to_average_fr.append(weight_WISE4*(params[13]/WISE4_synthetic_flux_erg))

		if type(fuv_galex_synthetic_flux_erg) != str:
			weight_fuv_galex = 1-(rel_err_fuv_galex/total_relative_err)
			to_average_fr.append(weight_fuv_galex*(params[16]/fuv_galex_synthetic_flux_erg))

		if type(nuv_galex_synthetic_flux_erg) != str:
			weight_nuv_galex = 1-(rel_err_nuv_galex/total_relative_err)
			to_average_fr.append(weight_nuv_galex*(params[17]/nuv_galex_synthetic_flux_erg))


		avg_flux_ratio = np.sum(to_average_fr)/len(to_average_fr)		#this is theta^2
		radius_in_parsec = np.sqrt(avg_flux_ratio*((params[4])**2))
		radius_in_solar_radii = radius_in_parsec/(2.256*10**(-8))


		to_chi_square = []
		absolute_error = []		#observed-calculated flux
		filters_used = []


		if type(fuv_galex_synthetic_flux_erg) != str:
			to_chi_square.append(weight_fuv_galex*(params[16]-(fuv_galex_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[16]-fuv_galex_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('fuv')

		if type(nuv_galex_synthetic_flux_erg) != str:
			to_chi_square.append(weight_nuv_galex*(params[17]-(nuv_galex_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[17]-nuv_galex_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('nuv')

		if type(sdss_u_synthetic_flux_erg) != str:
			to_chi_square.append(weight_sdss_u*(params[5]-(sdss_u_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[5]-sdss_u_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('u')

		if type(Johnson_B_synthetic_flux_erg) != str:
			to_chi_square.append(weight_Johnson_B*(params[14]-(Johnson_B_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[14]-Johnson_B_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('B')

		if type(sdss_g_synthetic_flux_erg) != str:
			to_chi_square.append(weight_sdss_g*(params[6]-(sdss_g_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[6]-sdss_g_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('g')

		if type(Johnson_V_synthetic_flux_erg) != str:
			to_chi_square.append(weight_Johnson_V*(params[15]-(Johnson_V_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[15]-Johnson_V_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('V')

		if type(sdss_r_synthetic_flux_erg) != str:
			to_chi_square.append(weight_sdss_r*(params[7]-(sdss_r_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[7]-sdss_r_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('r')

		if type(sdss_i_synthetic_flux_erg) != str:
			to_chi_square.append(weight_sdss_i*(params[8]-(sdss_i_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[8]-sdss_i_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('i')

		if type(sdss_z_synthetic_flux_erg) != str:
			to_chi_square.append(weight_sdss_z*(params[9]-(sdss_z_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[9]-sdss_z_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('z')

		if type(J_2MASS_synthetic_flux_erg) != str:
			to_chi_square.append(weight_J_2MASS*(params[18]-(J_2MASS_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[18]-J_2MASS_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('J')

		if type(H_2MASS_synthetic_flux_erg) != str:
			to_chi_square.append(weight_H_2MASS*(params[19]-(H_2MASS_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[19]-H_2MASS_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('H')

		if type(K_2MASS_synthetic_flux_erg) != str:
			to_chi_square.append(weight_K_2MASS*(params[20]-(K_2MASS_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[20]-K_2MASS_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('K')

		if type(WISE1_synthetic_flux_erg) != str:
			to_chi_square.append(weight_WISE1*(params[10]-(WISE1_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[10]-WISE1_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('W1')

		if type(WISE2_synthetic_flux_erg) != str:
			to_chi_square.append(weight_WISE2*(params[11]-(WISE2_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[11]-WISE2_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('W2')

		if type(WISE3_synthetic_flux_erg) != str:
			to_chi_square.append(weight_WISE3*(params[12]-(WISE3_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[12]-WISE3_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('W3')

		if type(WISE4_synthetic_flux_erg) != str:
			to_chi_square.append(weight_WISE4*(params[13]-(WISE4_synthetic_flux_erg*avg_flux_ratio))**2)
			absolute_error.append(params[13]-WISE4_synthetic_flux_erg*avg_flux_ratio)
			filters_used.append('W4')


		plt.plot(filters_used,absolute_error,'o')
		chi_square = np.sum(to_chi_square)

		avg_flux_ratio_plus_cs = (np.sum(to_average_fr)/len(to_average_fr))+chi_square
		radius_in_parsec_plus_cs = np.sqrt(avg_flux_ratio_plus_cs*((params[4])**2))
		radius_in_solar_radii_plus_cs = radius_in_parsec_plus_cs/(2.256*10**(-8))
		radius_error_in_solar_radii = radius_in_solar_radii_plus_cs-radius_in_solar_radii

		#mass = ((10**(LOGg))*((radius_in_solar_radii*(2.256*10**(-8)))**2)*((3.086*10**16)**2))/(1.9891*10**(30)*100*G)

		MASS_parameters_then_data = [chi_square,params[1],radius_in_solar_radii,params[0],params[2],params[3],params[4],radius_error_in_solar_radii,filters_used,absolute_error]
		return(MASS_parameters_then_data)


		#Temperature_values.append(Teff)
		#Temperature_err_values.append(rand_err_for_Teff)
		#Gravity_values.append(LOGg)
		#Gravity_err_values.append(rand_err_for_LOGg)
		#EBV_values.append(EBV)
		#EBV_err_values.append(rand_err_for_EBV)
		#Rv_values.append(Rv)
		#Rv_err_values.append(rand_err_for_Rv)
		#Distance_Pc_values.append(DISTANCE)
		#Distance_Pc_err_values.append(rand_err_for_DISTANCE)
		#Radius_values.append(radius_in_solar_radii)
		#Radius_err_values.append(radius_error_in_solar_radii)
		#chi_square_values.append(chi_square)
		#Mass_values.append(mass)
		#Mass_err_values.append()



		#print("temp: ", len(Temperature_values))
		#print("temp_err: ", len(Temperature_err_values))
		#print("grav: ", len(Gravity_values))
		#print("grav_err: ", len(Gravity_err_values))
		#print("ebv: ", len(EBV_values))
		#print("ebv_err: ", len(EBV_err_values))
		#print("Rv: ", len(Rv_values))
		#print("Rv_err: ", len(Rv_err_values))
		#print("distance: ", len(Distance_Pc_values))
		#print("distance_err: ", len(Distance_Pc_err_values))
		#print("radius: ", len(Radius_values))
		#print("radius_err: ", len(Radius_err_values))
		#print("chi square: ", len(chi_square_values))
		#print("mass: ", len(Mass_values))

if __name__ == '__main__':
        with Pool(7) as p:
                results = p.map(f,parameters)

chi_square_error = []
sed_radii = []
sed_logg = []

for list in results:
	chi_square_error.append(list[0])
	sed_radii.append(list[2])
	sed_logg.append(list[1])

#DO WEIGHTED CHI SQAURE



for i in range(len(chi_square_error)):
	if chi_square_error[i] == np.min(chi_square_error):
		mass = ((10**(sed_logg[i]))*((sed_radii[i]*(2.256*10**(-8)))**2)*((3.086*10**16)**2))/(1.9891*10**(30)*100*G)
		break

print('radius:',sed_radii[i])
print('mass:',mass)

os.chdir('/media/intern2/DD82-B92F/RESEARCH_FALL_2022/saved_sed_values/')

with open(star_name+'.txt', "w") as f:		#this method saves all sed parameters and results. It will name each file based on star_name, and will override data if a star_name.txt file already exists
	for i in range(0,len(results)):
		f.write(str(results[i][3]) + '\t' + str(results[i][1]) + '\t' + str(results[i][4]) + '\t' + str(results[i][5]) + '\t' + str(results[i][6]) + '\t' + str(results[i][2]) + '\t' + str(results[i][7]) + '\t' + str(results[i][0]) + '\t' + str(results[i][9]) + '\n')

with open(star_name+'_observed_params.txt', 'w') as f:
	f.write(str(temp_to_save) + '\t' + str(temp_err_to_save) + '\t' + str(grav_to_save) + '\t' + str(grav_err_to_save) + '\t' + str(ebv_to_save) + '\t' + str(ebv_err_to_save) + '\t' + str(rv_to_save) + '\t' + str(rv_err_to_save) + '\t' + str(distance_to_save) + '\t' + str(distance_err_to_save) + '\n')

					#make the absolute error vs filter plot
for l in range(len(sed_radii)):
	plt.plot(results[0][8],results[l][9],'o')

plt.title('Absolute Error vs. Filter: '+star_name)
plt.xlabel('Filter Used')
plt.ylabel('Observed Photometric flux - flux ratioed synthetic photometry flux')
plt.show()



				#make the radius histogram
bins = np.linspace(np.min(sed_radii),np.max(sed_radii),100)
plt.title('Radii Histogram: '+star_name)
plt.ylabel('# of Occurances')
plt.xlabel('Radii (solar radii)')
plt.hist(sed_radii,bins=bins)
plt.show()

#print('lowest chi_square radius in solar radii: ', Radius_values[final_count])
#print('lowest chi_square radius error in solar radii: ',Radius_err_values[final_count])
#print('lowest chi_square mass in solar masses: ', Mass_values[final_count])
#print('obs - calc err (last run only): ',absolute_error)
#print('lowest chi_square temp in K: ', Temperature_values[final_count])

#plt.plot(Radius_values,chi_square_values,'o')
#plt.title('Chi_squared vs Radius: '+star_name)
#plt.xlabel('Radius (in $\odot$)')
#plt.ylabel('Chi_squared value')
#plt.show()



