from transientclumps.TCGaussclumpsFunctions import *
from transientclumps.TCOffsetFunctions import *
from transientclumps.TCPrepFunctions import *
from transientclumps.merge_catalog import *
import glob as glob
import matplotlib.pyplot as plt

def pointing_check(direc,wave):

    sdf_files = sorted(glob.glob(direc+'/*'+wave+'*mJybmsm.sdf'))

    if wave == '450':
        pix_scale = 2.0
    elif wave == '850':
        pix_scale = 3.0

    region = sdf_files[0].split('/')[-1].split('_')[0]

    xoffs = []
    yoffs = []
    peakrs = []
    xoffuncs = []
    yoffuncs = []
    peakruncs = []
    datescans = []

    for eachfile in sdf_files:

        eachfile = eachfile

        datescans.append(eachfile.split('_000')[0].split('_')[-1]+'_000'+eachfile.split('_000')[-1].split('_')[0])

        prepare_image(eachfile)

        eachfile = eachfile.split('/')[-1]

        run_gaussclumps(eachfile.split('_mJybmsm')[0]+'_Jypbmsm_crop.sdf','parameters/GCParms.txt')

        target_catalogue = eachfile.split('_mJybmsm')[0]+'_Jypbmsm_crop'+'_log.FIT'
        reference_catalogue = sdf_files[0].split('/')[-1].split('_mJybmsm')[0]+'_Jypbmsm_crop'+'_log.FIT'

        offsets, errs, NSurviveCull, Nmatch, MatchInd, TargIndMatched, CulledInd, TargIndCulled  = source_match(target_catalogue, reference_catalogue, minpeak=0.2, maxrad=30, maxsep=10, cutoff=4, pix_scale=pix_scale)

        avg_xoff  = offsets[0]
        avg_yoff  = offsets[1]
        avg_peakr = offsets[2]
        std_xoff  = errs[0]
        std_yoff  = errs[1]
        std_peakr = errs[2]

        xoffs.append(avg_xoff)
        yoffs.append(avg_yoff)
        peakrs.append(avg_peakr)
        xoffuncs.append(std_xoff)
        yoffuncs.append(std_yoff)
        peakruncs.append(std_peakr)


        cat_match_name = target_catalogue.split('.FIT')[0]+'_match.FIT'

        merge_catalog(target_catalogue, reference_catalogue, eachfile.split('_000')[0].split('_')[-1],eachfile.split('_000')[-1].split('_')[0], cat_match_name,ref_index=MatchInd,cat_index=TargIndMatched)

        cat_cull_name = target_catalogue.split('.FIT')[0]+'_cull.FIT'
        
        merge_catalog(target_catalogue, reference_catalogue, eachfile.split('_000')[0].split('_')[-1],eachfile.split('_000')[-1].split('_')[0], cat_cull_name,ref_index=CulledInd,cat_index=TargIndCulled)




    os.system('mv *_Jypbmsm_crop* '+direc)
    plt.errorbar(xoffs,yoffs,xerr=xoffuncs,yerr=yoffuncs,linestyle='None',marker='o')
    plt.xlabel('RA Offsets (arcsecs)')
    plt.ylabel('Dec Offsets (arcsecs)')
    plt.savefig('pointsource_results/'+region+'/'+region+'_'+wave+'_local_pointing_corrections.pdf',format='pdf')

    plt.clf()

    plt.hist(peakrs)
    plt.xlabel('Flux Ratio with Reference')
    plt.ylabel('Count')
    plt.savefig('pointsource_results/'+region+'/'+region+'_'+wave+'_local_flux_ratio_with_ref.pdf',format='pdf')

    plt.clf()
