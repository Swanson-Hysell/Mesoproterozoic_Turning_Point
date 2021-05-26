# standard modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.cm import get_cmap
import matplotlib.patches as patches
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats

# pmagpy
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import cartopy.crs as ccrs
import cartopy
import xml.etree.ElementTree as ET
import pygplates as pgp

def create_vgp_FeatureCollection(compilation):
    """
    Loop through poles and produce a pygplates FeatureCollection of multiple VGPs.
    
    Modified from code by Michael G. Tetley.
    
    Parameters
    ----------
    compilation : dataframe
        pole compilation
        
    Returns
    -------
    vpgFeatureCollection : FeatureCollection
        pygplates FeatureCollection of VGPs in compilation
    """
    poleLat = []
    poleLon = []
    poleName = []
    poleSiteLat = []
    poleSiteLon = []
    poleNominalAge = []
    poleA95 = []
    poleAgeLowerLimit = []
    poleAgeUpperLimit = []
    plateID = []

    count = 0

    for i in range(len(compilation)):

        if np.isfinite(compilation['slat'][i]) and np.isfinite(compilation['slon'][i]) and \
           np.isfinite(compilation['age lower'][i]) and np.isfinite(compilation['age upper'][i]):

            poleLat.append(compilation['plat'][i])
            poleLon.append(compilation['plon'][i])

            poleName.append(compilation['name'][i] + ' (' + compilation['grade'][i] + ')')
            poleSiteLat.append(compilation['slat'][i])
            poleSiteLon.append(compilation['slon'][i])
            poleNominalAge.append(compilation['age'][i])
            poleA95.append(compilation['a95'][i])

            poleAgeLowerLimit.append(compilation['age lower'][i])
            poleAgeUpperLimit.append(compilation['age upper'][i])

            plateID.append(compilation['plateID'][i])

            count = count + 1

        # Print if any of the isfinite tests fail
        else:

            print('Bad data for : {}'.format(compilation['name'][i]))


    # Create new GPlates Feature Collection
    vpgFeatureCollection = pgp.FeatureCollection()

    # Create new GPlates feature 'VirtualGeomagneticPole'.
    # Pole lat, pole lon, pole name, and reconstruction plate ID added within PointOnSphere method.
    # Inc, Dec, A95, Age and Sample site lat/lon values to added within 'other_properties' method.

    for j in range(count):

        vgpFeature = pgp.Feature.create_reconstructable_feature(
                     pgp.FeatureType.create_gpml('VirtualGeomagneticPole'),
                     pgp.PointOnSphere([np.float(poleLat[j]), np.float(poleLon[j])]),
                     name = poleName[j],
                     reconstruction_plate_id = int(plateID[j]),
                     other_properties = [(pgp.PropertyName.create_gpml('poleA95'),
                                          pgp.XsDouble(np.float64(poleA95[j]))),
                                         (pgp.PropertyName.create_gpml('averageAge'),
                                          pgp.XsDouble(np.float64(poleNominalAge[j]))),
                                         (pgp.PropertyName.create_gpml('averageSampleSitePosition'),
                                          pgp.GmlPoint(pgp.PointOnSphere([np.float(poleSiteLat[j]), 
                                                                          np.float(poleSiteLon[j])])))])

        # Add newly created feature to existing Feature Collection
        vpgFeatureCollection.add(vgpFeature)

    return vpgFeatureCollection
    
def create_vgp_gpml(vpgFeatureCollection, filename):
    """
    Create a .gpml for a FeatureCollection of VGPs.
    
    Modified from code by Michael G. Tetley.
    
    Parameters
    ----------
    vpgFeatureCollection : FeatureCollection
        pygplates FeatureCollection of VGPs in compilation
        
    filename : string
        path and name for output .gpml
    """
    # Generate GPML output file
    gpmlOutputFile = filename

    # Check for existing output file with same name and remove if found
    if os.path.isfile(filename):
        os.remove(filename)

    # Check to make sure vgpFeatureCollection (feature collection) is not empty before writing to file
    if len(vpgFeatureCollection) != 0:
        outputFeatureCollection = pgp.FeatureCollectionFileFormatRegistry()
        outputFeatureCollection.write(vpgFeatureCollection, filename)

    # Check if new file was created and confirm export
    if os.path.isfile(filename):
        print('Palaeomagnetic pole data successfully exported in GPML format.')
        
def get_craton_XYs(gpml, plateIDs):
    """
    Get XY coordinates of a plate polygon from a .gpml.
    
    Parameters
    ----------
    gpml : string
        Path to .gpml file.
        
    plateIDs : list
        Of plateIDs.
    """
    # namespace dictionary
    ns = {'gpml':'http://www.gplates.org/gplates',
          'gml':'http://www.opengis.net/gml'}
    
    # initial parse
    tree = ET.parse(gpml)
    root = tree.getroot()
    
    # storage
    Xs = []
    Ys = []
    
    # iterate through featureMembers
    for featureMember in root.findall('gml:featureMember',ns):
        
        # get child
        for child in featureMember:
            slice_ind = child.tag.find('}')
            child_root = 'gpml:' + child.tag[slice_ind+1:]
        
        # check plateID
        plateID_path = child_root + '/gpml:reconstructionPlateId/gpml:ConstantValue/gpml:value'
        feature_plateID = int(featureMember.find(plateID_path,ns).text)
        if feature_plateID in plateIDs:
            
            if featureMember.find(child_root + '/gpml:outlineOf', ns)!=None:
                polygon_root = child_root + '/gpml:outlineOf'
            elif featureMember.find(child_root + '/gpml:boundary', ns)!=None:
                polygon_root = child_root + '/gpml:boundary'
            elif featureMember.find(child_root + '/gpml:unclassifiedGeometry', ns)!=None:
                polygon_root = child_root + '/gpml:unclassifiedGeometry'
            elif featureMember.find(child_root + '/gpml:centerLineOf', ns)!=None:
                polygon_root = child_root + '/gpml:centerLineOf'
            else:
                raise Exception('polygon_root undefined.')
            
            # get coordinates
            posList_path = polygon_root + '/gpml:ConstantValue/gpml:value/gml:Polygon/gml:exterior/gml:LinearRing/gml:posList'
            for feature_posList in featureMember.findall(posList_path,ns):
                np_posList = np.fromstring(feature_posList.text, dtype=float, sep=' ')
            
                # split into lat and lon
                lat_inds = np.arange(0, len(np_posList), 2, dtype=int)
                lon_inds = np.arange(1, len(np_posList), 2, dtype=int)

                feature_lat = np_posList[lat_inds]
                feature_lon = np_posList[lon_inds]
            
                Xs.append(feature_lon)
                Ys.append(feature_lat)
            
    return Xs, Ys

def get_single_craton_XYs(gpml):
    """
    Get XY coordinates of a plate polygon from a .gpml.
    
    Parameters
    ----------
    gpml : string
        Path to .gpml file.
    """
    # namespace dictionary
    ns = {'gpml':'http://www.gplates.org/gplates',
          'gml':'http://www.opengis.net/gml'}
    
    # initial parse
    tree = ET.parse(gpml)
    root = tree.getroot()
    
    # storage
    Xs = []
    Ys = []
    
    # iterate through featureMembers
    featureMember = root.find('gml:featureMember',ns)
        
    # get child
    for child in featureMember:
        slice_ind = child.tag.find('}')
        child_root = 'gpml:' + child.tag[slice_ind+1:]

    if featureMember.find(child_root + '/gpml:outlineOf', ns)!=None:
        polygon_root = child_root + '/gpml:outlineOf'
    elif featureMember.find(child_root + '/gpml:boundary', ns)!=None:
        polygon_root = child_root + '/gpml:boundary'
    elif featureMember.find(child_root + '/gpml:unclassifiedGeometry', ns)!=None:
        polygon_root = child_root + '/gpml:unclassifiedGeometry'
    elif featureMember.find(child_root + '/gpml:centerLineOf', ns)!=None:
        polygon_root = child_root + '/gpml:centerLineOf'
    else:
        raise Exception('polygon_root undefined.')

    # get coordinates
    posList_path = polygon_root + '/gpml:ConstantValue/gpml:value/gml:Polygon/gml:exterior/gml:LinearRing/gml:posList'
    for feature_posList in featureMember.findall(posList_path,ns):
        np_posList = np.fromstring(feature_posList.text, dtype=float, sep=' ')

        # split into lat and lon
        lat_inds = np.arange(0, len(np_posList), 2, dtype=int)
        lon_inds = np.arange(1, len(np_posList), 2, dtype=int)

        feature_lat = np_posList[lat_inds]
        feature_lon = np_posList[lon_inds]

        Xs = feature_lon
        Ys = feature_lat
            
    return Xs, Ys

def craton_plot(ax, plateIDs, Eulers, edgecolor, facecolor, alpha, linewidth, gpml = '../GPlates/Cratons/shapes_cratons.gpml', reverse_draw=False,draw_face=True,draw_edge=True):
    """
    Plot cratons with rotation.
    
    Parameters
    ----------
    ax : map axis
        On which to plot.
    
    plateIDs : list
        Of plateIDs.
        
    Eulers : list of lists
        Of Euler rotation parameters - if more than one given,
        the rotations will be additive.
    """
    # get cratons from .gpml
    
    Xs, Ys = get_craton_XYs(gpml, plateIDs)
    
    # draw in reverse
    if reverse_draw:
        Xs = np.flip(Xs)
        Ys = np.flip(Ys)
    
    # rotate cratons
    rotated_Xs = []
    rotated_Ys = []
    for i in range(len(Xs)):
        rotated_X = np.array([])
        rotated_Y = np.array([])
        for j in range(len(Xs[i])):
            this_X = [Xs[i][j]]
            this_Y = [Ys[i][j]]
            for k in range(len(Eulers)):
                this_Y, this_X = pmag.pt_rot(Eulers[k], this_Y, this_X)
            rotated_X = np.append(rotated_X, this_X)
            rotated_Y = np.append(rotated_Y, this_Y)
        rotated_Xs.append(rotated_X)
        rotated_Ys.append(rotated_Y)
        
    # add cratons
    for i in range(len(rotated_Xs)):
        XY = np.stack([rotated_Xs[i][::-1],rotated_Ys[i][::-1]],axis=1)
        if draw_edge:
            poly_edge = patches.Polygon(XY,
                                        edgecolor=edgecolor,facecolor='none',alpha=alpha,
                                        transform=ccrs.Geodetic(), linewidth=linewidth)
            ax.add_patch(poly_edge)
        if draw_face:
            poly_face = patches.Polygon(XY,
                                        edgecolor='none',facecolor=facecolor,alpha=alpha,
                                        transform=ccrs.Geodetic())
            ax.add_patch(poly_face)

        
        
def single_craton_plot(ax, gpml, Eulers, edgecolor, facecolor, alpha, linewidth):
    """
    Plot cratons with rotation.
    
    Parameters
    ----------
    ax : map axis
        On which to plot.
    
    gpml : string
        Path to .gpml file.
        
    Eulers : list of lists
        Of Euler rotation parameters - if more than one given,
        the rotations will be additive.
    """
    # get cratons from .gpml
    Xs, Ys = get_single_craton_XYs(gpml)
    
    # rotate craton
    rotated_Xs = np.array([])
    rotated_Ys = np.array([])
    for i in range(len(Xs)):
        this_X = [Xs[i]]
        this_Y = [Ys[i]]
        for j in range(len(Eulers)):
            this_Y, this_X = pmag.pt_rot(Eulers[j], this_Y, this_X)
        rotated_Xs = np.append(rotated_Xs, this_X)
        rotated_Ys = np.append(rotated_Ys, this_Y)
        
    # add craton
    #XY = np.stack([rotated_Xs[::-1],rotated_Ys[::-1]],axis=1)
    XY = np.stack([rotated_Xs,rotated_Ys],axis=1)
    print(XY.shape)
    poly_edge = patches.Polygon(XY,
                                edgecolor=edgecolor,facecolor='none',alpha=alpha,
                                transform=ccrs.Geodetic(), linewidth=linewidth)
    poly_face = patches.Polygon(XY,
                                edgecolor='none',facecolor=facecolor,alpha=alpha,
                                transform=ccrs.Geodetic())
    ax.add_patch(poly_face)
    ax.add_patch(poly_edge)
    
    
def equi_filled(map_axis, centerlon, centerlat, radius, color, alpha=1.0, edge_alpha=1.0):
    """
    Modified from the ipmag function equi().
    """
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = ipmag.shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])
    
    X = X[::-1]
    Y = Y[::-1]
    
    XY = np.stack([X,Y],axis=1)
    
    circle_edge = patches.Polygon(XY,
                                  edgecolor=color,facecolor='none',alpha=edge_alpha,
                                  transform=ccrs.Geodetic())
    circle_face = patches.Polygon(XY,
                                  edgecolor='none',facecolor=color,alpha=alpha,
                                  transform=ccrs.Geodetic())
    
    map_axis.add_patch(circle_face)
    map_axis.add_patch(circle_edge)
    
    
def rotated_pole_plot(ax, plon, plat, a95, Eulers, marker, s, marker_color, a95_color, a95_alpha, a95_edge_alpha=1.0):
    """
    Plot paleomagnetic pole with rotation.
    """
    # rotate pole
    rotated_plat = plat
    rotated_plon = plon
    for i in range(len(Eulers)):
        rotated_plat, rotated_plon = pmag.pt_rot(Eulers[i], [rotated_plat], [rotated_plon])
        rotated_plat = rotated_plat[0]
        rotated_plon = rotated_plon[0]
    
    # degrees to km conversion
    a95_km = a95 * 111.32
    
    # pole
    ax.scatter(rotated_plon, rotated_plat, marker=marker,
               color=marker_color, edgecolors='k', s=s,
               label='__nolegend__', zorder=101, transform=ccrs.PlateCarree())
    
    # a95
    equi_filled(ax, rotated_plon, rotated_plat, a95_km, a95_color, alpha=a95_alpha, edge_alpha=a95_edge_alpha)
    

def rotated_point_plot(ax, plon, plat, Eulers, s, marker_color):
    """
    Plot point with rotation.
    """
    # rotate pole
    rotated_plat = plat
    rotated_plon = plon
    for i in range(len(Eulers)):
        rotated_plat, rotated_plon = pmag.pt_rot(Eulers[i], rotated_plat, rotated_plon)
        rotated_plat = rotated_plat
        rotated_plon = rotated_plon
    
    # pole
    ax.scatter(rotated_plon, rotated_plat,
               color=marker_color, s=s, edgecolor='none',
               label='__nolegend__', zorder=101, transform=ccrs.PlateCarree())
    
def APWP_plot(ax,
              plon_rot, plat_rot, a95_rot, age_rot,
              plon_fix, plat_fix, a95_fix, age_fix,
              Euler,
              label_rot, color_rot, marker_rot, s_rot,
              label_fix, color_fix, marker_fix, s_fix,
              age_lim=None, cmap='viridis'):
    """
    Plot apparent polar wander paths for two connected cratons.
    """
    # get vmin and vmax, and slice data if necessary
    if age_lim==None:
        vmin = np.min([age_rot.min(), age_fix.min()])
        vmax = np.max([age_rot.max(), age_fix.max()])
    else:
        vmin = age_lim[0]
        vmax = age_lim[1]
        
        rot_mask = (age_rot>=age_lim[0])&(age_rot<=age_lim[1])
        plon_rot = plon_rot[rot_mask]
        plat_rot = plat_rot[rot_mask]
        a95_rot = a95_rot[rot_mask]
        age_rot = age_rot[rot_mask]
        
        fix_mask = (age_fix>=age_lim[0])&(age_fix<=age_lim[1])
        plon_fix = plon_fix[fix_mask]
        plat_fix = plat_fix[fix_mask]
        a95_fix = a95_fix[fix_mask]
        age_fix = age_fix[fix_mask]
    
    # rotate poles
    rotated_plat = np.array([])
    rotated_plon = np.array([])
    
    for i in range(len(plon_rot)):
        this_plat, this_plon = pmag.pt_rot(Euler, [plat_rot[i]], [plon_rot[i]])
        rotated_plat = np.append(rotated_plat, this_plat[0])
        rotated_plon = np.append(rotated_plon, this_plon[0])
    
    # degrees to km conversion
    a95_rot_km = a95_rot * 111.32
    a95_fix_km = a95_fix * 111.32
    
    # colormap to age
    color_mapping = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    colors_rot = color_mapping.to_rgba(age_rot).tolist()
    colors_fix = color_mapping.to_rgba(age_fix).tolist()
    
    # pole
    ax.scatter(plon_fix, plat_fix, marker=marker_fix,
               color=colors_fix, edgecolors='k', s=s_fix,
               label='__nolegend__', zorder=101, transform=ccrs.Geodetic())
    ax.scatter(rotated_plon, rotated_plat, marker=marker_rot,
               color=colors_rot, edgecolors='k', s=s_rot,
               label='__nolegend__', zorder=101, transform=ccrs.Geodetic())
    
    # a95
    for i in range(len(plon_fix)):
        equi_filled(ax, plon_fix[i], plat_fix[i], a95_fix_km[i], color_fix, alpha=0.3)
    for i in range(len(rotated_plon)):
        equi_filled(ax, rotated_plon[i], rotated_plat[i], a95_rot_km[i], color_rot, alpha=0.3)
        
    # create fake legend
    ax.scatter([], [], marker=marker_fix,
               color=color_fix, edgecolors='k', s=s_fix,
               label=label_fix, transform=ccrs.Geodetic())
    ax.scatter([], [], marker=marker_rot,
               color=color_rot, edgecolors='k', s=s_rot,
               label=label_rot, transform=ccrs.Geodetic())
    
    # colorbar
    color_mapping._A = []
    plt.colorbar(color_mapping, orientation='horizontal', shrink=0.8,
                 pad=0.05, label='age [Ma]')
    
    # prettify
    ax.legend(loc=4, fontsize=12)
    ax.set_title('{} rotated to {} - Euler : {}'.format(label_rot,label_fix,str(Euler)))
    
    plt.show()