# Author: Xiaoxuan Li
import arcpy
import pandas as pd
import os
import arcpy
from arcpy.sa import *
from datetime import date
import shutil
arcpy.env.overwriteOutput = True
import subprocess



def map_cover_agb_75m_finebeam():
    equation = r"E:\ScanSAR\ScanSAR\Result_0625\site_75m/equation.csv"
    df = pd.read_csv(equation)

    output_cover_dir = r"E:\ScanSAR\ScanSAR\Result_0906\ALS_predictions\FineBeam_cover/"
    output_agbd_dir =  r"E:\ScanSAR\ScanSAR\Result_0906\ALS_predictions\FineBeam_biomass/"
    chm = r"E:\ScanSAR\ScanSAR\Result_0906\SAR_FineBeam/FineBeam_20170702_75m.tif"

    suf = '20170706'
    print(suf)

    a_cover = df.loc[df['Date'] == int(suf), 'a'].values[0]
    b_cover = df.loc[df['Date'] == int(suf), 'b'].values[0]

    a_AGBD = df.loc[df['Date'] == int(suf), 'a'].values[1]
    b_AGBD = df.loc[df['Date'] == int(suf), 'b'].values[1]

    ras_cover = Exp(Raster(chm)*a_cover + b_cover)
    ras_cover.save(output_cover_dir + "FineBeam_" + suf + "_cover.tif")

    ras_AGBD = Exp(Raster(chm)*a_AGBD + b_AGBD)
    ras_AGBD.save(output_agbd_dir + "FineBeam_" + suf + "_agb.tif")


def extract_to_point_ScanSAR_GEDI():
    arc_dir = r"D:\temp\ArcGIS_project\ScanSAR_New\ScanSAR_New.gdb/"
    inpolygon = r"E:\ScanSAR\ScanSAR\Result_0906\GEDI_footprint/" + "GEDI_SAS.shp"
    ras_dir = r"E:\ScanSAR\ScanSAR\Result_0906\SAR/"
    chms = ['2014-09-06','2016-09-03','2017-09-02','2018-09-29','2019-09-28','2020-09-26','2021-09-25','2022-09-24']
    for chm in chms:
        print("Input ras: " + chm)
        chm = ras_dir + "PALSAR2_ScanSAR_HV_mtf_5_db_" + chm + ".tif"
        table = arc_dir + "ScanSAR_GEDI_" + \
                chm.split(".")[0].split("_")[7].split("-")[0] + \
                chm.split(".")[0].split("_")[7].split("-")[1] + \
                chm.split(".")[0].split("_")[7].split("-")[2]

        print("Output shp..." + table)
        ExtractValuesToPoints(inpolygon, chm, table, "INTERPOLATE", "VALUE_ONLY")


def table_ScanSAR_xlsx_GEDI():
    dir = r"D:\temp\ArcGIS_project\ScanSAR_New\ScanSAR_New.gdb/"
    out = r"E:\ScanSAR\ScanSAR\Result_0906\GEDI_footprint\SAR_Table/"
    shps = ['2014-09-06','2016-09-03','2017-09-02','2018-09-29','2019-09-28','2020-09-26','2021-09-25','2022-09-24']
    for shp in shps:
        print("Input ras: " + shp)
        ras = r"E:\ScanSAR\ScanSAR\single\PALSAR2_ScanSAR_HV_mtf_5_db_" + shp + ".tif"
        table = "ScanSAR_GEDI_" + \
                ras.split(".")[0].split("_")[6].split("-")[0] + \
                ras.split(".")[0].split("_")[6].split("-")[1] + \
                ras.split(".")[0].split("_")[6].split("-")[2]
        print("Output shp..." + table)
        arcpy.TableToTable_conversion(dir + table, out, table + ".csv")


def finebeam_GEDI():
    arc_dir = r"D:\temp\ArcGIS_project\ScanSAR_New\ScanSAR_New.gdb/"
    inpolygon = r"E:\ScanSAR\ScanSAR\Result_0906\GEDI_footprint/" + "GEDI_SAS.shp"
    out = r"E:\ScanSAR\ScanSAR\Result_0906\GEDI_footprint\SAR_Table"
    arcpy.env.workspace = r"E:\ScanSAR\ScanSAR\single/"
    chm = r"E:\ScanSAR\ScanSAR\Result_0906\SAR/FineBeam_20170702.tif"
    print("Input ras: " + chm)
    table = arc_dir + "FineBeam_GEDI_20170702"
    print("Output zonal stats..." + table)
    ExtractValuesToPoints(inpolygon, chm, table, "INTERPOLATE", "VALUE_ONLY")
    arcpy.TableToTable_conversion(table, out, "FineBeam_GEDI_20170702.csv")



def map_Cover_AGB_ScanSAR_L4A():
    equation = r"E:\ScanSAR\ScanSAR\Result_0906/GEDI_equation.csv"
    df = pd.read_csv(equation)
    output_agbd_dir =  r"E:\ScanSAR\ScanSAR\Result_0906\GEDI_predictions\ScanSAR_biomass/"
    arcpy.env.workspace = r"E:\ScanSAR\ScanSAR\Result_0906\SAR/"
    chms = arcpy.ListRasters('*.tif*')
    for chm in chms:
        if "PALSAR2_ScanSAR_HV" in chm:
            suf = chm.split("_")[6].split("-")[0] + chm.split("_")[6].split("-")[1] + chm.split("_")[6].split("-")[2].split(".")[0]
            if "20170706" != suf:
                print(suf)
                a_AGBD = df.loc[df['Date'] == int(suf), 'a'].values[0]
                b_AGBD = df.loc[df['Date'] == int(suf), 'b'].values[0]
                ras_AGBD = Exp(Raster(chm)*a_AGBD + b_AGBD)
                ras_AGBD.save(output_agbd_dir + "ScanSAR_GEDI_" + suf + "_agb.tif")


def map_Cover_AGB_finebeam_L4A():
    equation = r"E:\ScanSAR\ScanSAR\Result_0906/GEDI_equation.csv"
    df = pd.read_csv(equation)
    output_agbd_dir =  r"E:\ScanSAR\ScanSAR\Result_0906\GEDI_predictions\FineBeam_biomass/"
    arcpy.env.workspace = r"E:\ScanSAR\ScanSAR\Result_0906\SAR/"
    suf = "20170706"
    chm = r"E:\ScanSAR\ScanSAR\Result_0906\SAR/FineBeam_20170702.tif"
    print(suf)
    a_AGBD = df.loc[df['Date'] == int(suf), 'a'].values[0]
    b_AGBD = df.loc[df['Date'] == int(suf), 'b'].values[0]
    ras_AGBD = Exp(Raster(chm)*a_AGBD + b_AGBD)
    ras_AGBD.save(output_agbd_dir + "FineBeam_GEDI_" + suf + "_agb.tif")



def map_Cover_AGB_ScanSAR_SAS():
    equation = r"E:\ScanSAR\ScanSAR\Result_0906/GEDI_equation.csv"
    df = pd.read_csv(equation)
    output_agbd_dir =  r"E:\ScanSAR\ScanSAR\Result_0906\GEDI_SAS_predictions\ScanSAR_biomass/"
    arcpy.env.workspace = r"E:\ScanSAR\ScanSAR\Result_0906\SAR/"
    chms = arcpy.ListRasters('*.tif*')
    for chm in chms:
        if "PALSAR2_ScanSAR_HV" in chm:
            suf = chm.split("_")[6].split("-")[0] + chm.split("_")[6].split("-")[1] + chm.split("_")[6].split("-")[2].split(".")[0]
            if "20170706" != suf:
                print(suf)
                a_AGBD = df.loc[df['Date'] == int(suf), 'a'].values[1]
                b_AGBD = df.loc[df['Date'] == int(suf), 'b'].values[1]
                ras_AGBD = Exp(Raster(chm)*a_AGBD + b_AGBD)
                ras_AGBD.save(output_agbd_dir + "ScanSAR_SAS_" + suf + "_agb.tif")


def map_Cover_AGB_finebeam_SAS():
    equation = r"E:\ScanSAR\ScanSAR\Result_0906/GEDI_equation.csv"
    df = pd.read_csv(equation)
    output_agbd_dir =  r"E:\ScanSAR\ScanSAR\Result_0906\GEDI_SAS_predictions\FineBeam_biomass/"
    arcpy.env.workspace = r"E:\ScanSAR\ScanSAR\Result_0906\SAR/"
    suf = "20170706"
    chm = r"E:\ScanSAR\ScanSAR\Result_0906\SAR/FineBeam_20170702.tif"
    print(suf)
    a_AGBD = df.loc[df['Date'] == int(suf), 'a'].values[1]
    b_AGBD = df.loc[df['Date'] == int(suf), 'b'].values[1]
    ras_AGBD = Exp(Raster(chm)*a_AGBD + b_AGBD)
    ras_AGBD.save(output_agbd_dir + "FineBeam_SAS_" + suf + "_agb.tif")






def val_agbd_zonal():
    arc_dir = r"D:\temp\ArcGIS_project\ScanSAR_New\ScanSAR_New.gdb/"
    studysite = ['Agincourt', 'Welverdiendt', 'Justicia', 'Ireagh', 'Limpopo1', 'Limpopo2', 'Limpopo3']
    polygon_dir = r"E:\ScanSAR\ScanSAR\Result_1115\shp/"
    chm_full = r"E:\ScanSAR\ScanSAR\Result_1004\MC_100\Field_ALS_SAR/SAR_AGBD_2018_1ha.tif"
    out_dir = r"E:\ScanSAR\ScanSAR\Result_1115\table/"
    for site in studysite:
        shp = polygon_dir + site + "_1ha.shp"
        table = arc_dir + site + "_val_FAS_AGBD"
        out_table = site + "_val_FAS_AGBD.csv"
        ZonalStatisticsAsTable(shp, "FID", chm_full, table, "", "MEAN", "")
        arcpy.TableToTable_conversion(table, out_dir, out_table)
        print("Output shp..." + out_table)


def val_agbd_GLM_zonal():
    arc_dir = r"D:\temp\ArcGIS_project\ScanSAR_New\ScanSAR_New.gdb/"
    studysite = ['Agincourt', 'Welverdiendt', 'Justicia', 'Ireagh', 'Limpopo1', 'Limpopo2', 'Limpopo3']
    polygon_dir = r"E:\ScanSAR\ScanSAR\Result_1115\shp/"
    chm_full = r"E:\ScanSAR\ScanSAR\Result_1004\MC_100\Field_ALS_GEDI_SAR_glm/SAR_AGBD_2018_1ha.tif"
    out_dir = r"E:\ScanSAR\ScanSAR\Result_1115\table/"
    for site in studysite:
        shp = polygon_dir + site + "_1ha.shp"
        table = arc_dir + site + "_val_FAS_GLM_AGBD"
        out_table = site + "_val_FAS_GLM_AGBD.csv"
        ZonalStatisticsAsTable(shp, "FID", chm_full, table, "", "MEAN", "")
        arcpy.TableToTable_conversion(table, out_dir, out_table)
        print("Output shp..." + out_table)


def val_agbd_RF_zonal():
    arc_dir = r"D:\temp\ArcGIS_project\ScanSAR_New\ScanSAR_New.gdb/"
    studysite = ['Agincourt', 'Welverdiendt', 'Justicia', 'Ireagh', 'Limpopo1', 'Limpopo2', 'Limpopo3']
    polygon_dir = r"E:\ScanSAR\ScanSAR\Result_1115\shp/"
    chm_full = r"E:\ScanSAR\ScanSAR\Result_1004\MC_100\Field_ALS_GEDI_SAR_rf/SAR_AGBD_2018_1ha.tif"
    out_dir = r"E:\ScanSAR\ScanSAR\Result_1115\table/"
    for site in studysite:
        shp = polygon_dir + site + "_1ha.shp"
        table = arc_dir + site + "_val_FAS_RF_AGBD"
        out_table = site + "_val_FAS_RF_AGBD.csv"
        ZonalStatisticsAsTable(shp, "FID", chm_full, table, "", "MEAN", "")
        arcpy.TableToTable_conversion(table, out_dir, out_table)
        print("Output shp..." + out_table)


def val_H_zonal():
    arc_dir = r"D:\temp\ArcGIS_project\ScanSAR_New\ScanSAR_New.gdb/"
    studysite = ['Agincourt', 'Welverdiendt', 'Justicia', 'Ireagh', 'Limpopo1', 'Limpopo2', 'Limpopo3']
    polygon_dir = r"E:\ScanSAR\ScanSAR\Result_1115\shp/"
    out_dir = r"E:\ScanSAR\ScanSAR\Result_1115\table/"
    for site in studysite:
        chm = r"E:\ALS_archive\LiDAR_CHM_Mosaic_05/" + site + "_05.tif"
        shp = polygon_dir + site + "_1ha.shp"
        table = arc_dir + site + "_val_H"
        out_table = site + "_val_H.csv"
        ZonalStatisticsAsTable(shp, "FID", chm, table, "", "MEAN", "")
        arcpy.TableToTable_conversion(table, out_dir, out_table)
        print("Output shp..." + out_table)


def val_CC_zonal():
    arc_dir = r"D:\temp\ArcGIS_project\ScanSAR_New\ScanSAR_New.gdb/"
    studysite = ['Agincourt', 'Welverdiendt', 'Justicia', 'Ireagh', 'Limpopo1', 'Limpopo2', 'Limpopo3']
    polygon_dir = r"E:\ScanSAR\ScanSAR\Result_1115\shp/"
    out_dir = r"E:\ScanSAR\ScanSAR\Result_1115\table/"
    for site in studysite:
        chm = r"E:\ALS_archive\LiDAR_CHM_CC15/" + site + ".tif"
        shp = polygon_dir + site + "_1ha.shp"
        table = arc_dir + site + "_val_CC"
        out_table = site + "_val_CC.csv"
        ZonalStatisticsAsTable(shp, "FID", chm, table, "", "SUM", "")
        arcpy.TableToTable_conversion(table, out_dir, out_table)
        print("Output shp..." + out_table)


def resample_100_multi_ras():
    #dir = r"D:\temp\ArcGIS_project\ScanSAR_New\ScanSAR_New.gdb/"
    input = r"E:\ScanSAR\ScanSAR\single/"
    out = r"E:\ScanSAR\ScanSAR\Result_1231\Ras_100/"
    suf1 = ['2014-09-06','2016-09-03','2017-09-02','2018-09-29','2019-09-28','2020-09-26','2021-09-25','2022-09-24']
    suf = ['2023-09-23']
    for ras in suf:
        print("Ras: " + ras)
        ras_input = input + "PALSAR2_ScanSAR_HV_mtf_5_db_" + ras + ".tif"
        ras_output = out + "Ras_" + ras.split("-")[0] + "_100.tif"
        arcpy.Resample_management(ras_input, ras_output, "100", "BILINEAR")


def resample_clip():
    shp_dir = r"D:\temp\ArcGIS_project\ScanSAR_New\ScanSAR_New.gdb\KNPboundary_Project"
    out_dir = r"E:\ScanSAR\ScanSAR\Result_1231\Ras_100_clip/"
    arcpy.MakeFeatureLayer_management(shp_dir, 'shp')
    desc = arcpy.Describe("shp")
    extent = str(desc.extent.XMin) + " " + \
             str(desc.extent.YMin) + " " + \
             str(desc.extent.XMax) + " " + \
             str(desc.extent.YMax)
    arcpy.env.workspace = r"E:\ScanSAR\ScanSAR\Result_1231\Ras_100/"
    chms = arcpy.ListRasters('*.tif*')
    for chm in chms:
        if "2023" in chm:
            print("Output..." + chm)
            arcpy.Clip_management(chm, extent, out_dir+chm,'shp', -999, "ClippingGeometry", "")


def resample_clip_zonal():
    arc_dir = r"D:\temp\ArcGIS_project\ScanSAR_New\ScanSAR_New.gdb/"
    inpolygon = arc_dir + "GEDI_KNP_SAR_100"
    ras_dir = r"E:\ScanSAR\ScanSAR\Result_1231\Ras_100_clip/"
    suf1 = ['2014-09-06','2016-09-03','2017-09-02','2018-09-29','2019-09-28','2020-09-26','2021-09-25','2022-09-24']
    suf = ['2023-09-23']
    for ras in suf:
        ras_SAR = ras_dir + "Ras_" + ras.split("-")[0] + "_100.tif"
        table = arc_dir + "KNP_100_" + ras.split("-")[0]
        print(inpolygon)
        print(ras_SAR)
        print(table)
        ExtractValuesToPoints(inpolygon, ras_SAR, table, "", "VALUE_ONLY")
        out_dir = r"E:\ScanSAR\ScanSAR\Result_1231\Table_zonal/"
        outTable = "KNP_SAR_100_" + ras.split("-")[0] + ".csv"
        arcpy.TableToTable_conversion(table, out_dir, outTable)
        print("Output shp..." + outTable)



from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from sys import argv

def MCH_CC_AGBD():
    dir = r"E:\Biomass\CSIR\Biomass_2018\CSIR_25m.csv"
    df = pd.read_csv(dir)

    X = np.array(df['MEAN05'])
    Y = np.array(df['CC15'])
    Z = np.array(df['agbd_ha'])
    fig = plt.figure()
    plt.figure(figsize=(10, 20), dpi=300)
    plt.rcParams.update({'font.size': 4})
    plt.scatter(X, Y, s=4, c=Z, cmap='viridis')
    plt.colorbar().ax.set_title("AGBD (Mg/ha)")
    plt.xlabel("MCH (m)")
    plt.ylabel("CC (%)")
    plt.show()


def MCH_AGBD():
    dir = r"E:\Biomass\CSIR\Biomass_2018\CSIR_25m.csv"
    df = pd.read_csv(dir)
    X = np.array(df['MEAN05'])
    Y = np.array(df['agbd_ha'])
    Z = np.array(df['CC15'])
    fig = plt.figure()
    plt.figure(figsize=(10, 10), dpi=300)
    plt.rcParams.update({'font.size': 4})
    plt.scatter(X, Y, s=4, c=Z, cmap='viridis')
    plt.colorbar().ax.set_title("CC (%)")
    plt.ylabel("AGBD (Mg/ha)")
    plt.xlabel("MCH (m)")
    plt.show()





def MCH_CC_histplot():
    import numpy as np
    import seaborn as sns
    from matplotlib import pyplot as plt
    import matplotlib.pylab as pylab
    fig, ax1 = plt.subplots()
    plt.rcParams.update({'font.size': 4})
    dir = r"E:\Biomass\CSIR\Biomass_2018\CSIR_25m.csv"
    df = pd.read_csv(dir)
    X_field = np.array(df['MEAN05'])
    Y_field = np.array(df['CC15'])
    dir = r"E:\Biomass\CSIR\Result_01292024\merge.csv"
    df = pd.read_csv(dir)
    X = np.array(df['MCH'])
    Y = np.array(df['CC'])
    sns.histplot(x=X, y=Y, cbar=True)
    ax1.scatter(X_field, Y_field, s=4, c="red")
    cbar = ax1.collections[0].colorbar
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label("Density", fontsize=10, labelpad=-40, y=1.05, rotation=0)
    ax1.set_xlabel("MCH (m)")
    ax1.set_ylabel("CC (%)")
    plt.show()

