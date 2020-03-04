# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is code from ~2016-2018.
# It is used to carry out an MCMC (emcee) in order to estimate the parameters of a biconical outflow.
# The input is the velocity and velocity error of the two components of a double-peaked profile.
# The output is the position of the walkers at various points in time.
# The parameters of the bicone are phi, the inclination, theta, the PA on the sky, r_t, the turnover radius
# for the velocity law, the half opening angle and the maximum velocity.
# ~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                                                                                                     
import matplotlib as plt
from pylab import *
import numpy as np
from decimal import *
import numpy
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib.pyplot as plt
#from astropy.modeling import models, fitting                                                                                                                                                      
from pylab import *
import math
from mpl_toolkits.mplot3d import Axes3D
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TIME TO DO THE BICONE                                                                                                                                            






# These fxns construct the various bicones based on the input parameters
# The output is the 3-dimensional construction of x, y, and z coordinates and the LOS velocity
# at that point.

# There are three different 'construction' fxns, one for the bicone,
# one for the cocone, and one for the nesting cone, all of which
# have different geometries, see Nevin et al. 2018 for more
# details about how I build these.


# This code constructs a three dimensional bicone that will be fed into lnlike_bicone 
# where it will be projected down to 2D
def bicone_construction(phi,theta,r_t, half_angle, vel_max):
    print('these are the input parameters',phi,theta,r_t, half_angle, vel_max)
    # the five input parameters determine the geometry of the bicone,
    # phi is the inclination, theta, is the position angle on the sky in 2D,
    # r_t is the turnover radius, or point at which the material moving along
    # the walls begins to decelerate,
    # half_angle, is the half opening angle,
    # and vel_max is the maximum intrinsic (meaning 3D) velocity along the walls
    # of the cone.
    
    # The model only moves material along the walls of the cone, so it could
    # be called an 'evacuated bicone'
    
    half_angle_r=math.radians(half_angle)
    phi_1=math.radians(phi)
    theta_1=math.radians(theta)
    psi_1=math.radians(0)

    # the total heigh is twice the effective radius because
    # the velocity decays linearly from vmax to zero
    h=2*r_t
    # to get the radius of the cone, you need to do some trig
    # this is the radius across the opening
    r=math.tan(half_angle_r)*h
    
    # this is your sampling around the cone, here I'm taking 50 points
    theta_para=(np.asarray(np.linspace(0,2*np.pi,30)))
    
    # This is the sampling along the heigh of the cone
    u1_orig=np.asarray(np.linspace(-h,0,10))+h#
    u2_orig=-np.asarray(np.linspace(-h,0,10))-h
    u=-np.asarray(np.linspace(-h,0,10))

    


    

    #R is the rotation matrix
    R=np.matrix([[np.cos(theta_1)*np.cos(psi_1),np.cos(phi_1)*np.sin(psi_1)+np.sin(phi_1)*np.sin(theta_1)*np.cos(psi_1),
                  np.sin(phi_1)*np.sin(psi_1)-np.cos(phi_1)*np.sin(theta_1)*np.cos(psi_1)],
                     [-np.cos(theta_1)*np.sin(psi_1),
                      np.cos(phi_1)*np.cos(psi_1)-np.sin(phi_1)*np.sin(theta_1)*np.sin(psi_1),
                      np.sin(phi_1)*np.cos(psi_1)+np.cos(phi_1)*np.sin(theta_1)*np.sin(psi_1)],
                     [np.sin(theta_1), -np.sin(phi_1)*np.cos(theta_1), np.cos(phi_1)*np.cos(theta_1)]])
    #x is going to be the front facing side (blueshifted) and x_2 is the rear cone
    x=[]
    x_2=[]
    y=[]
    y_2=[]
    z=[]
    z_2=[]
    
    vel_radial_bottom=[]
    vel_radial_top=[]
    
    x_totes=[]
    y_totes=[]
    z_totes=[]
    x_plane=[]
    y_plane=[]
    z_plane=[]
    x_plane_2=[]
    y_plane_2=[]
    z_plane_2=[]
    
    from mpl_toolkits import mplot3d
    plt.clf()
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(0,0,0,marker='*', color='orange')
    
    for i in range(len(theta_para)):
        for j in range(len(u)):
            if i==0:
                i=1
            if j==0:
                j=1

            #here's the parametric equation for a bicone
            X=((h-u[j])*r*np.cos((theta_para[i])))/h
            Y=((h-u[j])*r*np.sin((theta_para[i])))/h

          

            #blue cone above is Z_1
            #red cone below is Z_2
            Z_1=u2_orig[j]#was u2
            Z_2=u1_orig[j]#was u1
            
            
            
            y_plane_into=0

            pos=np.matrix([X,Y,Z_1])
            pos_2=np.matrix([X,Y,Z_2])
            
            # so the 'plane' positions are all at y=0, because they are projected
            # onto the plane of the sky
            pos_plane=np.matrix([X,y_plane_into,Z_1])
            pos_plane_2=np.matrix([X,y_plane_into, Z_2])
            
            # take the dot product of the rotation matrix and the current position
            # in order to move to the location on the cone
            new_pos=np.dot(R,pos.transpose())
            new_pos_2=np.dot(R,pos_2.transpose())
            new_pos_plane=np.dot(R,pos_plane.transpose())
            x_plane.append(new_pos_plane[0])
            
            y_plane.append(new_pos_plane[1])
            z_plane.append(new_pos_plane[2])
            new_pos_plane_2=np.dot(R, pos_plane_2.transpose())
            x_plane_2.append(new_pos_plane_2[0])
            y_plane_2.append(new_pos_plane_2[1])
            z_plane_2.append(new_pos_plane_2[2])
            
            
            #Z_2=u1[j]
            x.append(new_pos[0])
            x_2.append(new_pos_2[0])
            y.append(new_pos[1])
            y_2.append(new_pos_2[1])
            
            z.append(new_pos[2])
            z_2.append(new_pos_2[2])
            
            
            if np.sqrt(np.array(new_pos[0]).ravel()**2+np.array(new_pos[1]).ravel()**2+np.array(new_pos[2]).ravel()**2)<r_t:
                # so if the point is at a distance less than the turnover radius, then the velocity is the max velocity
                vel_max_top=-vel_max
            if np.sqrt(np.array(new_pos[0]).ravel()**2+np.array(new_pos[1]).ravel()**2+np.array(new_pos[2]).ravel()**2)>r_t:
                # if the point is greater than the turnover radius, then it will decline linearly to zero at 2*r_t
                vel_max_top=-(vel_max-(vel_max/r_t)*(np.sqrt(np.array(new_pos[0]).ravel()**2+np.array(new_pos[1]).ravel()**2+np.array(new_pos[2]).ravel()**2)-r_t))
                
            if np.sqrt(np.array(new_pos[0]).ravel()**2+np.array(new_pos[1]).ravel()**2+np.array(new_pos[2]).ravel()**2)>2*r_t:
                # if the point is farther than 2*r_t, then it is already at zero
                vel_max_top=0
            #do that thing:
            #equation: d/dw/e((x^2 + y^2)/c^2 = (z-z_0)^2)
            #recall c=r/h
            #<2*x,2*y,-(2*r^2*z)/h^2>
            #plug in everything and normalize and then the z component is your v_rad
    
            # Time to project this velocity into a line of sight velocity
            slope_x=np.array(new_pos[0]).ravel()#eliminated 2*
            slope_y=np.array(new_pos[1]).ravel()# eliminated 2*
            slope_z=np.array(new_pos[2]).ravel()#-(2*r**2*np.array(new_pos[2]).ravel())/h**2
            
            norm=np.sqrt(slope_x**2+slope_y**2+slope_z**2)
            A_comp=-(slope_y/norm)
            ax.plot3D([slope_x,0],[slope_y,0],[slope_z,0])
            
            
            # multiply the intrinsic velocity by the component of A that is in the y direction
            # y, which is the line of sight direction
            # currently, the velocity is projected along the wall of the bicone
            rad_proj=vel_max_top*np.array(A_comp).ravel()
            
            
            #print('position', new_pos)
            #print('total length to x,y,z = 0', norm, 'y cont', slope_y)
            #print('velocity at this point', vel_max_top)
            #print('component', *np.array(A_comp).ravel())
            #print('projected vel', rad_proj)
            
            vel_radial_top.append(round(rad_proj[0]))
            
            slope_x=2*np.array(new_pos_2[0]).ravel()
            slope_y=2*np.array(new_pos_2[1]).ravel()
            
            slope_z_2=-(2*r**2*np.array(new_pos_2[2]).ravel())/h**2
            norm=np.sqrt(slope_x**2+slope_y**2+slope_z_2**2)
            
            if np.sqrt(np.array(new_pos_2[0]).ravel()**2+np.array(new_pos_2[1]).ravel()**2+np.array(new_pos_2[2]).ravel()**2)<r_t:
                
            
#                vel_max_bottom=-(vel_max/r_t)*np.sqrt(new_pos_2[0]**2+new_pos_2[1]**2+new_pos_2[2]**2)
                vel_max_bottom=-vel_max
            if np.sqrt(np.array(new_pos_2[0]).ravel()**2+np.array(new_pos_2[1]).ravel()**2+np.array(new_pos_2[2]).ravel()**2)>r_t:
                #vel_max_bottom=vel_max-(vel_max/r_t)*np.sqrt(new_pos_2[0]**2+new_pos_2[1]**2+new_pos_2[2]**2)
                vel_max_bottom=-(vel_max-(vel_max/r_t)*(np.sqrt(np.array(new_pos[0]).ravel()**2+np.array(new_pos[1]).ravel()**2+np.array(new_pos[2]).ravel()**2)-r_t))
                #if new_pos_2[2]>0:
                #    vel_max_bottom=-vel_max_bottom
                #if vel_max_bottom>0:
                #    vel_max_bottom=0
            if np.sqrt(np.array(new_pos_2[0]).ravel()**2+np.array(new_pos_2[1]).ravel()**2+np.array(new_pos_2[2]).ravel()**2)>2*r_t:
                #vel_max_bottom=vel_max-(vel_max/r_t)*np.sqrt(new_pos_2[0]**2+new_pos_2[1]**2+new_pos_2[2]**2)
                vel_max_bottom=0
            #if abs(vel_max_bottom)>1001:
            #    vel_max_bottom=0
            #A dot x should give us the vector projection along x - divide by the magnitude of A to get cos (theta)
            A_comp=-(slope_y/norm)
            
            
            
            rad_proj_2=np.array(A_comp).ravel()*vel_max_bottom
            vel_radial_bottom.append(round(rad_proj_2[0]))
    
    # planes at x=0, y=0, and z=0 to guide the eye
    xx, zz = np.meshgrid(range(-20,20), range(-20,20))
    im = ax.scatter3D(x, y, z, c=vel_radial_top, cmap='Blues_r')
    im2 = ax.scatter3D(x_2,y_2,z_2,c=vel_radial_bottom, cmap='Reds')
    # calculate corresponding z
    yy = [0 for x in xx]
    ax.plot_surface(xx, yy, zz, alpha=0.2)
    #ax.set_ylim([-20,0])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.view_init(azim=0)
    plt.colorbar(im)
    plt.colorbar(im2)
    plt.savefig('3D_diag.png')
    
    
    # So the first half of all these merged arrays is the blueshifted side of the cone (above the plane of the sky)
    # and x_2, etc are the redshifted, or 'bottom' side of the cone.
    merged_x=x+x_2
    merged_y=y+y_2
    merged_z=z+z_2
    
    merged_x_plane=x_plane+x_plane_2
    merged_y_plane=y_plane+y_plane_2
    merged_z_plane=z_plane+z_plane_2
    
    
    
    merged_vel=vel_radial_top+vel_radial_bottom
    
    
    plt.clf()
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    im = ax.scatter3D(merged_x_plane, merged_y_plane, merged_z_plane, c=merged_vel, cmap='RdBu_r')
    plt.colorbar(im)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.plot_surface(xx, yy, zz, alpha=0.2)
    ax.view_init(azim=10)
    
    plt.savefig('3D_merged_plane.png')
    
    plt.clf()
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    im = ax.scatter3D(merged_x, merged_y, merged_z, c=merged_vel, cmap='RdBu_r')
    plt.colorbar(im)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.plot_surface(xx, yy, zz, alpha=0.2)
    ax.view_init(azim=10)
    plt.savefig('3D_merged.png')
    STOP
    
    return merged_x, merged_y, merged_z, merged_vel, merged_x_plane, merged_y_plane, merged_z_plane
                
    
            
# This creates a diagnostic plot by returning the 
def diagnostic_bicone_plot(z, xs_data_concat, xs_data_concat_ortho):


    #bicone_construction(h,phi,theta,r_t, half_angle )
    phi_array, theta_array, r_t,half_angle, vel_max_array=z


    h=2*r_t
    #h=[2*x for x in r_t]

    out=bicone_construction(phi_array, theta_array, r_t, half_angle, vel_max_array)

    # So this gives you a whole lot of 3D measurements that you do not actually need

    merged_x=[float(y) for y in out[0]]
    merged_y=[float(y) for y in out[1]]
    merged_z=[float(y) for y in out[2]]
    merged_vel=[float(y) for y in out[3]]
    merged_x_plane=[float(y) for y in out[4]]
    merged_y_plane=[float(y) for y in out[5]]
    merged_z_plane=[float(y) for y in out[6]]

    '''cone_red_x=[float(y) for y in out[7]]
    cone_red_y=[float(y) for y in out[8]]
    cone_red_z=[float(y) for y in out[9]]
    cone_blue_x=[float(y) for y in out[10]]
    cone_blue_y=[float(y) for y in out[11]]
    cone_blue_z=[float(y) for y in out[12]]

    cone_red_vel=[float(y) for y in out[13]]
    cone_blue_vel=[float(y) for y in out[14]]'''
    ##vel_radial_top actually refers to the stuff going away from us on the bottom -ugh
    #it's the projection looking down in z
    #p=ax.scatter(x_2,y_2,z_2, c=vel_radial_bottom)


    #%```````````````````````````````````
    #I'm going to make the phi cuts according to observe PAs
    PA_obs_1=np.radians(PA1+90)
    PA_obs_2=np.radians(PA2+90)

    '''x_slit_blue=np.linspace(-100*h,100*h, len(cone_blue_z))
    x_slit_red=np.linspace(-100*h,100*h,len(cone_red_z))'''
    
    x_slit=np.linspace(-100*h,100*h, len(merged_x))

    len_pts=max(upper_lim-lower_lim,upper_lim_ortho-lower_lim_ortho)

    #height is based upon the length of the overall data we have on the bicone and the slitwidth
    height=slitwidth/pixelscale_1#(slitwidth/(pixelscale*len_pts))#*(h)#I kinda want this to be a fraction of something else#10#was 4

    if PA1 < 90:
        z_slit_1_upper=[x*np.tan(PA_obs_1)-height/np.cos(PA_obs_1) for x in merged_x]
        z_slit_1_lower=[x*np.tan(PA_obs_1)+height/np.cos(PA_obs_1) for x in merged_x]
        
    else:
        z_slit_1_upper=[x*np.tan(PA_obs_1)+height/np.cos(PA_obs_1) for x in merged_x]
        z_slit_1_lower=[x*np.tan(PA_obs_1)-height/np.cos(PA_obs_1) for x in merged_x]
        
    height=slitwidth/pixelscale_2  
    if PA2 > 90:
        z_slit_2_upper=[x*np.tan(PA_obs_2)+height/np.cos(PA_obs_2) for x in merged_x]
        z_slit_2_lower=[x*np.tan(PA_obs_2)-height/np.cos(PA_obs_2) for x in merged_x]
        
    else: 
        z_slit_2_upper=[x*np.tan(PA_obs_2)-height/np.cos(PA_obs_2) for x in merged_x]
        z_slit_2_lower=[x*np.tan(PA_obs_2)+height/np.cos(PA_obs_2) for x in merged_x]
        
   

    if PA1<90:
        inds_1=np.logical_and(np.logical_and((np.asarray(merged_z) > np.asarray(z_slit_1_lower)), (np.asarray(merged_z) < np.asarray(z_slit_1_upper))), (np.asarray(merged_y) < np.asarray(merged_y_plane)))
        inds_1_below=np.logical_and(np.logical_and((np.asarray(merged_z) > np.asarray(z_slit_1_lower)), (np.asarray(merged_z) < np.asarray(z_slit_1_upper))), (np.asarray(merged_y) > np.asarray(merged_y_plane)))
    else:
        inds_1=np.logical_and(np.logical_and((np.asarray(merged_z) < np.asarray(z_slit_1_lower)), (np.asarray(merged_z) > np.asarray(z_slit_1_upper))), (np.asarray(merged_y) < np.asarray(merged_y_plane)))
        inds_1_below=np.logical_and(np.logical_and((np.asarray(merged_z) < np.asarray(z_slit_1_lower)), (np.asarray(merged_z) > np.asarray(z_slit_1_upper))), (np.asarray(merged_y) > np.asarray(merged_y_plane)))
    if PA2>90:
        inds_2=np.logical_and(np.logical_and((np.asarray(merged_z) < np.asarray(z_slit_2_lower)), (np.asarray(merged_z) > np.asarray(z_slit_2_upper))), (np.asarray(merged_y) < np.asarray(merged_y_plane)))
        inds_2_below=np.logical_and(np.logical_and((np.asarray(merged_z) < np.asarray(z_slit_2_lower)), (np.asarray(merged_z) > np.asarray(z_slit_2_upper))), (np.asarray(merged_y) > np.asarray(merged_y_plane)))
    else:
        inds_2=np.logical_and(np.logical_and((np.asarray(merged_z) > np.asarray(z_slit_2_lower)), (np.asarray(merged_z) < np.asarray(z_slit_2_upper))), (np.asarray(merged_y) < np.asarray(merged_y_plane)))
        inds_2_below=np.logical_and(np.logical_and((np.asarray(merged_z) > np.asarray(z_slit_2_lower)), (np.asarray(merged_z) < np.asarray(z_slit_2_upper))), (np.asarray(merged_y) > np.asarray(merged_y_plane)))
    

    '''xs_cont_1 is the blueshifted cone of PA1'''
    xs_cont_1=np.asarray(merged_x)[np.asarray(inds_1)]
    ys_cont_1=np.asarray(merged_y)[np.asarray(inds_1)]
    zs_cont_1=np.asarray(merged_z)[np.asarray(inds_1)]
    vel_cont_1=np.asarray(merged_vel)[np.asarray(inds_1)]

    '''xs_cont_1_below is the redshifted cone of PA1'''
    xs_cont_1_below=np.asarray(merged_x)[np.asarray(inds_1_below)]
    ys_cont_1_below=np.asarray(merged_y)[np.asarray(inds_1_below)]
    zs_cont_1_below=np.asarray(merged_z)[np.asarray(inds_1_below)]
    vel_cont_1_below=np.asarray(merged_vel)[np.asarray(inds_1_below)]

    xs_cont_2=np.asarray(merged_x)[np.asarray(inds_2)]
    ys_cont_2=np.asarray(merged_y)[np.asarray(inds_2)]
    zs_cont_2=np.asarray(merged_z)[np.asarray(inds_2)]
    vel_cont_2=np.asarray(merged_vel)[np.asarray(inds_2)]
    xs_cont_2_below=np.asarray(merged_x)[np.asarray(inds_2_below)]
    ys_cont_2_below=np.asarray(merged_y)[np.asarray(inds_2_below)]
    zs_cont_2_below=np.asarray(merged_z)[np.asarray(inds_2_below)]
    vel_cont_2_below=np.asarray(merged_vel)[np.asarray(inds_2_below)]

    
   
    
    r_slit_1_above=[]
    v_slit_1_above=[]
    r_slit_2_above=[]
    v_slit_2_above=[]
    r_slit_1_below=[]
    v_slit_1_below=[]
    r_slit_2_below=[]
    v_slit_2_below=[]

    mask_slit_1_above=[]
    mask_slit_2_above=[]
    mask_slit_1_below=[]
    mask_slit_2_below=[]


    stepping=1
    
    
    



    length=max((upper_lim-lower_lim),(upper_lim_ortho-lower_lim_ortho))

    if PA1 < 90:
        additive=np.pi/2
    else:
        additive=np.pi/2
    for j in range(length):
        #
        width=pixelscale_1*j
        '''z_1_perp_upper=[-x*(1/np.tan(PA_obs_1))-(width+stepping+0.5)/np.sin(math.radians(PA1-90)) for x in xs_cont_1]
        z_1_perp_lower=[-x*(1/np.tan(PA_obs_1))-(width+0.5)/np.sin(math.radians(PA1-90)) for x in xs_cont_1]
        z_1_perpl_upper=[-x*(1/np.tan(PA_obs_1))+(width+0.5)/np.sin(math.radians(PA1-90)) for x in xs_cont_1]
        z_1_perpl_lower=[-x*(1/np.tan(PA_obs_1))+(width+stepping+0.5)/np.sin(math.radians(PA1-90)) for x in xs_cont_1]
        '''
        z_1_perp_upper=[-x*(1/np.tan(PA_obs_1))+(width+1.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1]
        z_1_perp_lower=[-x*(1/np.tan(PA_obs_1))+(width+0.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1]
        z_1_perpl_upper=[-x*(1/np.tan(PA_obs_1))-(width+0.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1]
        z_1_perpl_lower=[-x*(1/np.tan(PA_obs_1))-(width+1.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1]


        if PA1 < 90:
            pos_1=np.logical_and((np.asarray(zs_cont_1) > np.asarray(z_1_perp_upper)), (np.asarray(zs_cont_1) < np.asarray(z_1_perp_lower)))
            pos_1_1=np.logical_and((np.asarray(zs_cont_1) > np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1) < np.asarray(z_1_perpl_lower)))
            pos_1_center=np.logical_and((np.asarray(zs_cont_1) < np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1) > np.asarray(z_1_perp_lower)))
        else:
            pos_1=np.logical_and((np.asarray(zs_cont_1) < np.asarray(z_1_perp_upper)), (np.asarray(zs_cont_1) > np.asarray(z_1_perp_lower)))
            pos_1_1=np.logical_and((np.asarray(zs_cont_1) < np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1) > np.asarray(z_1_perpl_lower)))
            pos_1_center=np.logical_and((np.asarray(zs_cont_1) > np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1) < np.asarray(z_1_perp_lower)))
        
        
       

        if j==0:
            xs_cont_slit_middle=np.asarray(xs_cont_1)[np.asarray(pos_1_center)]
            zs_cont_slit_middle=np.asarray(zs_cont_1)[np.asarray(pos_1_center)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit_middle)**2+np.mean(zs_cont_slit_middle)**2))==True:
                r_slit_1_above.append(j)
                v_slit_1_above.append(10)
                mask_slit_1_above.append(1)

            else:
                #so this is the positions to the North so we'll make them negative
                r_slit_1_above.append(j)
                v_slit_1_above.append(np.mean(np.asarray(vel_cont_1)[np.asarray(pos_1_center)]))
                mask_slit_1_above.append(0)
        else:
            xs_cont_slit=np.asarray(xs_cont_1)[np.asarray(pos_1)]
            zs_cont_slit=np.asarray(zs_cont_1)[np.asarray(pos_1)]


            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))==True:

                r_slit_1_above.append(j)
                v_slit_1_above.append(10)
                mask_slit_1_above.append(1)
            else:
                #so this is the positions to the North so we'll make them negative
                r_slit_1_above.append(j)
                v_slit_1_above.append(np.mean(np.asarray(vel_cont_1)[np.asarray(pos_1)]))
                mask_slit_1_above.append(0)
                #r_slit_1_above.append(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))
                #v_slit_1_above.append(np.mean(np.asarray(vel_cont_1)[np.asarray(pos_1)]))


            xs_cont_slit=np.asarray(xs_cont_1)[np.asarray(pos_1_1)]
            zs_cont_slit=np.asarray(zs_cont_1)[np.asarray(pos_1_1)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))==True:

                r_slit_1_above.append(-j)
                v_slit_1_above.append(10)
                mask_slit_1_above.append(1)
            else:
                r_slit_1_above.append(-j)
                #r_slit_1_above.append(-np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))
                #r_slit_1_above.append(-np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))
                v_slit_1_above.append(np.mean(np.asarray(vel_cont_1)[np.asarray(pos_1_1)]))
                mask_slit_1_above.append(0)

        z_1_perp_upper=[-x*(1/np.tan(PA_obs_1))+(width+1.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1_below]
        z_1_perp_lower=[-x*(1/np.tan(PA_obs_1))+(width+0.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1_below]
        z_1_perpl_upper=[-x*(1/np.tan(PA_obs_1))-(width+0.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1_below]
        z_1_perpl_lower=[-x*(1/np.tan(PA_obs_1))-(width+1.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1_below]
        
        if PA1 < 90:
            pos_1_below=np.logical_and((np.asarray(zs_cont_1_below) > np.asarray(z_1_perp_upper)), (np.asarray(zs_cont_1_below) < np.asarray(z_1_perp_lower)))
            pos_1_below_1=np.logical_and((np.asarray(zs_cont_1_below) > np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1_below) < np.asarray(z_1_perpl_lower)))
            pos_1_below_center=np.logical_and((np.asarray(zs_cont_1_below) < np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1_below) > np.asarray(z_1_perp_lower)))
        else:
            pos_1_below=np.logical_and((np.asarray(zs_cont_1_below) < np.asarray(z_1_perp_upper)), (np.asarray(zs_cont_1_below) > np.asarray(z_1_perp_lower)))
            pos_1_below_1=np.logical_and((np.asarray(zs_cont_1_below) < np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1_below) > np.asarray(z_1_perpl_lower)))
            pos_1_below_center=np.logical_and((np.asarray(zs_cont_1_below) > np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1_below) < np.asarray(z_1_perp_lower)))

        if j==0:
            xs_cont_slit_middle=np.asarray(xs_cont_1_below)[np.asarray(pos_1_below_center)]
            zs_cont_slit_middle=np.asarray(zs_cont_1_below)[np.asarray(pos_1_below_center)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit_middle)**2+np.mean(zs_cont_slit_middle)**2))==True:
                r_slit_1_below.append(j)
                v_slit_1_below.append(10)
                mask_slit_1_below.append(1)

            else:
                #so this is the positions to the North so we'll make them negative
                r_slit_1_below.append(j)
                v_slit_1_below.append(np.mean(np.asarray(vel_cont_1_below)[np.asarray(pos_1_below_center)]))
                mask_slit_1_below.append(0)
        else:

            xs_cont_slit=np.asarray(xs_cont_1_below)[np.asarray(pos_1_below)]
            zs_cont_slit=np.asarray(zs_cont_1_below)[np.asarray(pos_1_below)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))==True:
                r_slit_1_below.append(j)
                v_slit_1_below.append(10)
                mask_slit_1_below.append(1)

            else:
                r_slit_1_below.append(j)
                v_slit_1_below.append(np.mean(np.asarray(vel_cont_1_below)[np.asarray(pos_1_below)]))
                mask_slit_1_below.append(0)

            xs_cont_slit=np.asarray(xs_cont_1_below)[np.asarray(pos_1_below_1)]
            zs_cont_slit=np.asarray(zs_cont_1_below)[np.asarray(pos_1_below_1)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2)):
                r_slit_1_below.append(-j)
                v_slit_1_below.append(10)
                mask_slit_1_below.append(1)

            else:
                r_slit_1_below.append(-j)
                v_slit_1_below.append(np.mean(np.asarray(vel_cont_1_below)[np.asarray(pos_1_below_1)]))
                mask_slit_1_below.append(0)

        '''Now PA 2 which is greater than 90'''
        width=pixelscale_2*j

        z_2_perp_upper=[-x*(1/np.tan(PA_obs_2))+(width+1.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2]
        z_2_perp_lower=[-x*(1/np.tan(PA_obs_2))+(width+0.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2]
        z_2_perpl_upper=[-x*(1/np.tan(PA_obs_2))-(width+0.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2]
        z_2_perpl_lower=[-x*(1/np.tan(PA_obs_2))-(width+1.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2]
        
        
        if PA2 > 90:
            pos_1_2=np.logical_and((np.asarray(zs_cont_2) < np.asarray(z_2_perp_upper)), (np.asarray(zs_cont_2) > np.asarray(z_2_perp_lower)))
            pos_1_3=np.logical_and((np.asarray(zs_cont_2) < np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2)> np.asarray(z_2_perpl_lower)))

            pos_1_2_center=np.logical_and((np.asarray(zs_cont_2) > np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2) < np.asarray(z_2_perp_lower)))
        else:
            pos_1_2=np.logical_and((np.asarray(zs_cont_2) > np.asarray(z_2_perp_upper)), (np.asarray(zs_cont_2) < np.asarray(z_2_perp_lower)))
            pos_1_3=np.logical_and((np.asarray(zs_cont_2) > np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2)< np.asarray(z_2_perpl_lower)))

            pos_1_2_center=np.logical_and((np.asarray(zs_cont_2) < np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2) > np.asarray(z_2_perp_lower)))
        #
        
        #if j==0:
            #plt.plot(xs_cont_2, z_2_perp_upper, label='z_1_perp_upper',lw=4, color='black')
            #plt.plot(xs_cont_2, z_2_perp_lower, label='z_2_perp_lower',lw=4,ls='--', color='red')
            #plt.plot(xs_cont_2, z_2_perpl_upper, label='z_2_perpl_upper',lw=4, color='pink')
            #plt.show()
        if j==0:
            xs_cont_slit_middle=np.asarray(xs_cont_2)[np.asarray(pos_1_2_center)]
            zs_cont_slit_middle=np.asarray(zs_cont_2)[np.asarray(pos_1_2_center)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit_middle)**2+np.mean(zs_cont_slit_middle)**2))==True:

                r_slit_2_above.append(j)
                v_slit_2_above.append(10)
                mask_slit_2_above.append(1)

            else:
                #so this is the positions to the North so we'll make them negative
                r_slit_2_above.append(j)
                v_slit_2_above.append(np.mean(np.asarray(vel_cont_2)[np.asarray(pos_1_2_center)]))
                mask_slit_2_above.append(0)
        else:

            xs_cont_slit=np.asarray(xs_cont_2)[np.asarray(pos_1_2)]
            zs_cont_slit=np.asarray(zs_cont_2)[np.asarray(pos_1_2)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))==True:
                r_slit_2_above.append(j)
                v_slit_2_above.append(10)
                mask_slit_2_above.append(1)

            else:
                r_slit_2_above.append(j)
                v_slit_2_above.append(np.mean(np.asarray(vel_cont_2)[np.asarray(pos_1_2)]))
                mask_slit_2_above.append(0)

            xs_cont_slit=np.asarray(xs_cont_2)[np.asarray(pos_1_3)]
            zs_cont_slit=np.asarray(zs_cont_2)[np.asarray(pos_1_3)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))==True:

                r_slit_2_above.append(-j)
                v_slit_2_above.append(10)
                mask_slit_2_above.append(1)

            else:
                r_slit_2_above.append(-j)
                v_slit_2_above.append(np.mean(np.asarray(vel_cont_2)[np.asarray(pos_1_3)]))
                mask_slit_2_above.append(0)
            #math.isnan(np.mean(np.asarray(vel_cont_2)[np.asarray(pos_1_3)])):



        z_2_perp_upper=[-x*(1/np.tan(PA_obs_2))+(width+1.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2_below]
        z_2_perp_lower=[-x*(1/np.tan(PA_obs_2))+(width+0.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2_below]
        z_2_perpl_upper=[-x*(1/np.tan(PA_obs_2))-(width+0.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2_below]
        z_2_perpl_lower=[-x*(1/np.tan(PA_obs_2))-(width+1.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2_below]

        
        if PA2 >  90:
            pos_1_4=np.logical_and((np.asarray(zs_cont_2_below) < np.asarray(z_2_perp_upper)), (np.asarray(zs_cont_2_below) > np.asarray(z_2_perp_lower)))
            pos_1_5=np.logical_and((np.asarray(zs_cont_2_below) < np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2_below) > np.asarray(z_2_perpl_lower)))
            pos_1_4_center=np.logical_and((np.asarray(zs_cont_2_below) > np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2_below) < np.asarray(z_2_perp_lower)))
        else:
            pos_1_4=np.logical_and((np.asarray(zs_cont_2_below) > np.asarray(z_2_perp_upper)), (np.asarray(zs_cont_2_below) < np.asarray(z_2_perp_lower)))
            pos_1_5=np.logical_and((np.asarray(zs_cont_2_below) > np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2_below) < np.asarray(z_2_perpl_lower)))
            pos_1_4_center=np.logical_and((np.asarray(zs_cont_2_below) < np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2_below) > np.asarray(z_2_perp_lower)))


        if j==0:
            xs_cont_slit_middle=np.asarray(xs_cont_2_below)[np.asarray(pos_1_4_center)]
            zs_cont_slit_middle=np.asarray(zs_cont_2_below)[np.asarray(pos_1_4_center)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit_middle)**2+np.mean(zs_cont_slit_middle)**2))==True:

                r_slit_2_below.append(j)
                v_slit_2_below.append(10)
                mask_slit_2_below.append(1)
            else:
                #so this is the positions to the North so we'll make them negative
                r_slit_2_below.append(j)
                v_slit_2_below.append(np.mean(np.asarray(vel_cont_2_below)[np.asarray(pos_1_4_center)]))
                mask_slit_2_below.append(0)
        else:

            xs_cont_slit=np.asarray(xs_cont_2_below)[np.asarray(pos_1_4)]
            zs_cont_slit=np.asarray(zs_cont_2_below)[np.asarray(pos_1_4)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))==True:
                r_slit_2_below.append(j)
                v_slit_2_below.append(10)
                mask_slit_2_below.append(1)

            else:
                r_slit_2_below.append(j)
                v_slit_2_below.append(np.mean(np.asarray(vel_cont_2_below)[np.asarray(pos_1_4)]))
                mask_slit_2_below.append(0)

            xs_cont_slit=np.asarray(xs_cont_2_below)[np.asarray(pos_1_5)]
            zs_cont_slit=np.asarray(zs_cont_2_below)[np.asarray(pos_1_5)]

            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))==True:

                r_slit_2_below.append(-j)
                v_slit_2_below.append(10)
                mask_slit_2_below.append(1)

            else:
                r_slit_2_below.append(-j)
                v_slit_2_below.append(np.mean(np.asarray(vel_cont_2_below)[np.asarray(pos_1_5)]))
                mask_slit_2_below.append(0)



    inds = np.array(r_slit_1_above).argsort()
    sortedx_1_above = np.array(r_slit_1_above)[inds]#-minus_one
    sortedv_1_above = np.array(v_slit_1_above)[inds]

    sorted_mask = np.array(mask_slit_1_above)[inds]

    m_sortedx_1_above = ma.masked_array(sortedx_1_above, mask=sorted_mask)
    m_sortedv_1_above = ma.masked_array(sortedv_1_above, mask=sorted_mask)

    inds = np.array(r_slit_1_below).argsort()
    sortedx_1_below = np.array(r_slit_1_below)[inds]
    sortedv_1_below = np.array(v_slit_1_below)[inds]

    sorted_mask = np.array(mask_slit_1_below)[inds]

    m_sortedx_1_below = ma.masked_array(sortedx_1_below, mask=sorted_mask)
    m_sortedv_1_below = ma.masked_array(sortedv_1_below, mask=sorted_mask)

    inds = np.array(r_slit_2_above).argsort()
    sortedx_2_above = np.array(r_slit_2_above)[inds]
    sortedv_2_above = np.array(v_slit_2_above)[inds]

    sorted_mask = np.array(mask_slit_2_above)[inds]

    m_sortedx_2_above = ma.masked_array(sortedx_2_above, mask=sorted_mask)
    m_sortedv_2_above = ma.masked_array(sortedv_2_above, mask=sorted_mask)

    inds = np.array(r_slit_2_below).argsort()
    sortedx_2_below = np.array(r_slit_2_below)[inds]
    sortedv_2_below = np.array(v_slit_2_below)[inds]

    sorted_mask = np.array(mask_slit_2_below)[inds]

    m_sortedx_2_below = ma.masked_array(sortedx_2_below, mask=sorted_mask)
    m_sortedv_2_below = ma.masked_array(sortedv_2_below, mask=sorted_mask)



    '''Okay so we need to cut this according to the xs_data_concat array'''
    first=xs_data[lower_lim:upper_lim][0]
    last=xs_data[lower_lim:upper_lim][-1]
    for z in range(len(m_sortedx_1_below)):
        if sortedx_1_below[z]==first:
            i_b_first=z
        if sortedx_1_below[z]==last:
            i_b_last=z+1

    for z in range(len(m_sortedx_1_above)):
        if sortedx_1_above[z]==first:
            i_a_first=z
        if sortedx_1_above[z]==last:
            i_a_last=z+1


    xs_fit_concat=np.ma.concatenate([m_sortedx_1_below[i_b_first:i_b_last],m_sortedx_1_above[i_a_first:i_a_last]])#xs_fit_bin+xs_fit_bin_below
    vs_fit_concat=np.ma.concatenate([m_sortedv_1_below[i_b_first:i_b_last],m_sortedv_1_above[i_a_first:i_a_last]])#vs_fit_bin+vs_fit_bin_below

    first=xs_data_ortho[lower_lim_ortho:upper_lim_ortho][0]
    last=xs_data_ortho[lower_lim_ortho:upper_lim_ortho][-1]

    '''find index where m_sortedx_1_below and same for above have the values of first and last and cut it there'''
    for z in range(len(m_sortedx_2_below)):
        if sortedx_2_below[z]==first:
            i_b_first=z
        if sortedx_2_below[z]==last:
            i_b_last=z+1

    for z in range(len(m_sortedx_2_above)):
        if sortedx_2_above[z]==first:
            i_a_first=z
        if sortedx_2_above[z]==last:
            i_a_last=z+1




    xs_fit_concat_ortho=np.ma.concatenate([m_sortedx_2_below[i_b_first:i_b_last],m_sortedx_2_above[i_a_first:i_a_last]])#xs_fit_ortho_bin+xs_fit_ortho_bin_below
    vs_fit_concat_ortho=np.ma.concatenate([m_sortedv_2_below[i_b_first:i_b_last],m_sortedv_2_above[i_a_first:i_a_last]])#vs_fit_ortho_bin+vs_fit_ortho_bin_below

    
    return xs_fit_concat, vs_fit_concat, xs_fit_concat_ortho, vs_fit_concat_ortho    
    
    

def lnlike_bicone_three_slits(z, PA1, PA2, PA3, pixelscale_1, slitwidth):


    #bicone_construction(h,phi,theta,r_t, half_angle )
    phi_array, theta_array, r_t,half_angle, vel_max_array=z


    h=2*r_t
    #h=[2*x for x in r_t]

    out=bicone_construction(phi_array, theta_array, r_t, half_angle, vel_max_array)

    # So this gives you a whole lot of 3D measurements that you do not actually need

    merged_x=[float(y) for y in out[0]]
    merged_y=[float(y) for y in out[1]]
    merged_z=[float(y) for y in out[2]]
    merged_vel=[float(y) for y in out[3]]
    merged_x_plane=[float(y) for y in out[4]]
    merged_y_plane=[float(y) for y in out[5]]
    merged_z_plane=[float(y) for y in out[6]]

    plt.clf()
    plt.scatter(merged_x, merged_y, merged_z, c=merged_vel)
    plt.colorbar()
    plt.savefig('3D_cone.png')



    plt.clf()
    plt.scatter(merged_x_plane, merged_y_plane, merged_z_plane, c=merged_vel)
    plt.colorbar()
    plt.savefig('3D_cone_plane.png')
    
    
    #%```````````````````````````````````
    #I'm going to make the phi cuts according to observe PAs
    PA_obs_1=np.radians(PA1+90)
    PA_obs_2=np.radians(PA2+90)
    PA_obs_3=np.radians(PA3+90)

    x_slit=np.linspace(-100*h,100*h, len(merged_x))

    len_pts=max(upper_lim-lower_lim,upper_lim_ortho-lower_lim_ortho)

    #height is based upon the length of the overall data we have on the bicone and the slitwidth
    height=slitwidth/pixelscale_1#

    if PA1 < 90:
        z_slit_1_upper=[x*np.tan(PA_obs_1)-height/np.cos(PA_obs_1) for x in merged_x]
        z_slit_1_lower=[x*np.tan(PA_obs_1)+height/np.cos(PA_obs_1) for x in merged_x]
        
    else:
        z_slit_1_upper=[x*np.tan(PA_obs_1)+height/np.cos(PA_obs_1) for x in merged_x]
        z_slit_1_lower=[x*np.tan(PA_obs_1)-height/np.cos(PA_obs_1) for x in merged_x]
        
    height=slitwidth/pixelscale_2  
    if PA2 > 90:
        z_slit_2_upper=[x*np.tan(PA_obs_2)+height/np.cos(PA_obs_2) for x in merged_x]
        z_slit_2_lower=[x*np.tan(PA_obs_2)-height/np.cos(PA_obs_2) for x in merged_x]
        
    else: 
        z_slit_2_upper=[x*np.tan(PA_obs_2)-height/np.cos(PA_obs_2) for x in merged_x]
        z_slit_2_lower=[x*np.tan(PA_obs_2)+height/np.cos(PA_obs_2) for x in merged_x]
        
   

    if PA1<90:
        inds_1=np.logical_and(np.logical_and((np.asarray(merged_z) > np.asarray(z_slit_1_lower)), (np.asarray(merged_z) < np.asarray(z_slit_1_upper))), (np.asarray(merged_y) < np.asarray(merged_y_plane)))
        inds_1_below=np.logical_and(np.logical_and((np.asarray(merged_z) > np.asarray(z_slit_1_lower)), (np.asarray(merged_z) < np.asarray(z_slit_1_upper))), (np.asarray(merged_y) > np.asarray(merged_y_plane)))
    else:
        inds_1=np.logical_and(np.logical_and((np.asarray(merged_z) < np.asarray(z_slit_1_lower)), (np.asarray(merged_z) > np.asarray(z_slit_1_upper))), (np.asarray(merged_y) < np.asarray(merged_y_plane)))
        inds_1_below=np.logical_and(np.logical_and((np.asarray(merged_z) < np.asarray(z_slit_1_lower)), (np.asarray(merged_z) > np.asarray(z_slit_1_upper))), (np.asarray(merged_y) > np.asarray(merged_y_plane)))
    if PA2>90:
        inds_2=np.logical_and(np.logical_and((np.asarray(merged_z) < np.asarray(z_slit_2_lower)), (np.asarray(merged_z) > np.asarray(z_slit_2_upper))), (np.asarray(merged_y) < np.asarray(merged_y_plane)))
        inds_2_below=np.logical_and(np.logical_and((np.asarray(merged_z) < np.asarray(z_slit_2_lower)), (np.asarray(merged_z) > np.asarray(z_slit_2_upper))), (np.asarray(merged_y) > np.asarray(merged_y_plane)))
    else:
        inds_2=np.logical_and(np.logical_and((np.asarray(merged_z) > np.asarray(z_slit_2_lower)), (np.asarray(merged_z) < np.asarray(z_slit_2_upper))), (np.asarray(merged_y) < np.asarray(merged_y_plane)))
        inds_2_below=np.logical_and(np.logical_and((np.asarray(merged_z) > np.asarray(z_slit_2_lower)), (np.asarray(merged_z) < np.asarray(z_slit_2_upper))), (np.asarray(merged_y) > np.asarray(merged_y_plane)))
    

    # xs_cont_1 is the blueshifted cone of PA1
    xs_cont_1=np.asarray(merged_x)[np.asarray(inds_1)]
    ys_cont_1=np.asarray(merged_y)[np.asarray(inds_1)]
    zs_cont_1=np.asarray(merged_z)[np.asarray(inds_1)]
    vel_cont_1=np.asarray(merged_vel)[np.asarray(inds_1)]

    # xs_cont_1_below is the redshifted cone of PA1
    xs_cont_1_below=np.asarray(merged_x)[np.asarray(inds_1_below)]
    ys_cont_1_below=np.asarray(merged_y)[np.asarray(inds_1_below)]
    zs_cont_1_below=np.asarray(merged_z)[np.asarray(inds_1_below)]
    vel_cont_1_below=np.asarray(merged_vel)[np.asarray(inds_1_below)]

    xs_cont_2=np.asarray(merged_x)[np.asarray(inds_2)]
    ys_cont_2=np.asarray(merged_y)[np.asarray(inds_2)]
    zs_cont_2=np.asarray(merged_z)[np.asarray(inds_2)]
    vel_cont_2=np.asarray(merged_vel)[np.asarray(inds_2)]
    xs_cont_2_below=np.asarray(merged_x)[np.asarray(inds_2_below)]
    ys_cont_2_below=np.asarray(merged_y)[np.asarray(inds_2_below)]
    zs_cont_2_below=np.asarray(merged_z)[np.asarray(inds_2_below)]
    vel_cont_2_below=np.asarray(merged_vel)[np.asarray(inds_2_below)]

    
    

    plt.clf()
    plt.scatter(xs_cont_1, ys_cont_1, c=vel_cont_1)
    plt.colorbar()
    plt.savefig('diag_slit1.png')


    plt.clf()
    plt.scatter(xs_cont_2, ys_cont_2, c=vel_cont_2)
    plt.colorbar()
    plt.savefig('diag_slit2.png')
    STOP

    
    r_slit_1_above=[]
    v_slit_1_above=[]
    r_slit_2_above=[]
    v_slit_2_above=[]
    r_slit_1_below=[]
    v_slit_1_below=[]
    r_slit_2_below=[]
    v_slit_2_below=[]

    mask_slit_1_above=[]
    mask_slit_2_above=[]
    mask_slit_1_below=[]
    mask_slit_2_below=[]


    stepping=1
    
    
    



    length=max((upper_lim-lower_lim),(upper_lim_ortho-lower_lim_ortho))

    if PA1 < 90:
        additive=np.pi/2
    else:
        additive=np.pi/2
    for j in range(length):
        #
        width=pixelscale_1*j
        
        z_1_perp_upper=[-x*(1/np.tan(PA_obs_1))+(width+1.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1]
        z_1_perp_lower=[-x*(1/np.tan(PA_obs_1))+(width+0.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1]
        z_1_perpl_upper=[-x*(1/np.tan(PA_obs_1))-(width+0.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1]
        z_1_perpl_lower=[-x*(1/np.tan(PA_obs_1))-(width+1.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1]


        if PA1 < 90:
            pos_1=np.logical_and((np.asarray(zs_cont_1) > np.asarray(z_1_perp_upper)), (np.asarray(zs_cont_1) < np.asarray(z_1_perp_lower)))
            pos_1_1=np.logical_and((np.asarray(zs_cont_1) > np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1) < np.asarray(z_1_perpl_lower)))
            pos_1_center=np.logical_and((np.asarray(zs_cont_1) < np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1) > np.asarray(z_1_perp_lower)))
        else:
            pos_1=np.logical_and((np.asarray(zs_cont_1) < np.asarray(z_1_perp_upper)), (np.asarray(zs_cont_1) > np.asarray(z_1_perp_lower)))
            pos_1_1=np.logical_and((np.asarray(zs_cont_1) < np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1) > np.asarray(z_1_perpl_lower)))
            pos_1_center=np.logical_and((np.asarray(zs_cont_1) > np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1) < np.asarray(z_1_perp_lower)))
        
        
       

        if j==0:
            xs_cont_slit_middle=np.asarray(xs_cont_1)[np.asarray(pos_1_center)]
            zs_cont_slit_middle=np.asarray(zs_cont_1)[np.asarray(pos_1_center)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit_middle)**2+np.mean(zs_cont_slit_middle)**2))==True:
                r_slit_1_above.append(j)
                v_slit_1_above.append(10)
                mask_slit_1_above.append(1)

            else:
                #so this is the positions to the North so we'll make them negative
                r_slit_1_above.append(j)
                v_slit_1_above.append(np.mean(np.asarray(vel_cont_1)[np.asarray(pos_1_center)]))
                mask_slit_1_above.append(0)
        else:
            xs_cont_slit=np.asarray(xs_cont_1)[np.asarray(pos_1)]
            zs_cont_slit=np.asarray(zs_cont_1)[np.asarray(pos_1)]


            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))==True:

                r_slit_1_above.append(j)
                v_slit_1_above.append(10)
                mask_slit_1_above.append(1)
            else:
                #so this is the positions to the North so we'll make them negative
                r_slit_1_above.append(j)
                v_slit_1_above.append(np.mean(np.asarray(vel_cont_1)[np.asarray(pos_1)]))
                mask_slit_1_above.append(0)
                #r_slit_1_above.append(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))
                #v_slit_1_above.append(np.mean(np.asarray(vel_cont_1)[np.asarray(pos_1)]))


            xs_cont_slit=np.asarray(xs_cont_1)[np.asarray(pos_1_1)]
            zs_cont_slit=np.asarray(zs_cont_1)[np.asarray(pos_1_1)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))==True:

                r_slit_1_above.append(-j)
                v_slit_1_above.append(10)
                mask_slit_1_above.append(1)
            else:
                r_slit_1_above.append(-j)
                #r_slit_1_above.append(-np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))
                #r_slit_1_above.append(-np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))
                v_slit_1_above.append(np.mean(np.asarray(vel_cont_1)[np.asarray(pos_1_1)]))
                mask_slit_1_above.append(0)

        z_1_perp_upper=[-x*(1/np.tan(PA_obs_1))+(width+1.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1_below]
        z_1_perp_lower=[-x*(1/np.tan(PA_obs_1))+(width+0.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1_below]
        z_1_perpl_upper=[-x*(1/np.tan(PA_obs_1))-(width+0.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1_below]
        z_1_perpl_lower=[-x*(1/np.tan(PA_obs_1))-(width+1.5*pixelscale_1)/np.cos(PA_obs_1+additive) for x in xs_cont_1_below]
        
        if PA1 < 90:
            pos_1_below=np.logical_and((np.asarray(zs_cont_1_below) > np.asarray(z_1_perp_upper)), (np.asarray(zs_cont_1_below) < np.asarray(z_1_perp_lower)))
            pos_1_below_1=np.logical_and((np.asarray(zs_cont_1_below) > np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1_below) < np.asarray(z_1_perpl_lower)))
            pos_1_below_center=np.logical_and((np.asarray(zs_cont_1_below) < np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1_below) > np.asarray(z_1_perp_lower)))
        else:
            pos_1_below=np.logical_and((np.asarray(zs_cont_1_below) < np.asarray(z_1_perp_upper)), (np.asarray(zs_cont_1_below) > np.asarray(z_1_perp_lower)))
            pos_1_below_1=np.logical_and((np.asarray(zs_cont_1_below) < np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1_below) > np.asarray(z_1_perpl_lower)))
            pos_1_below_center=np.logical_and((np.asarray(zs_cont_1_below) > np.asarray(z_1_perpl_upper)), (np.asarray(zs_cont_1_below) < np.asarray(z_1_perp_lower)))

        if j==0:
            xs_cont_slit_middle=np.asarray(xs_cont_1_below)[np.asarray(pos_1_below_center)]
            zs_cont_slit_middle=np.asarray(zs_cont_1_below)[np.asarray(pos_1_below_center)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit_middle)**2+np.mean(zs_cont_slit_middle)**2))==True:
                r_slit_1_below.append(j)
                v_slit_1_below.append(10)
                mask_slit_1_below.append(1)

            else:
                #so this is the positions to the North so we'll make them negative
                r_slit_1_below.append(j)
                v_slit_1_below.append(np.mean(np.asarray(vel_cont_1_below)[np.asarray(pos_1_below_center)]))
                mask_slit_1_below.append(0)
        else:

            xs_cont_slit=np.asarray(xs_cont_1_below)[np.asarray(pos_1_below)]
            zs_cont_slit=np.asarray(zs_cont_1_below)[np.asarray(pos_1_below)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))==True:
                r_slit_1_below.append(j)
                v_slit_1_below.append(10)
                mask_slit_1_below.append(1)

            else:
                r_slit_1_below.append(j)
                v_slit_1_below.append(np.mean(np.asarray(vel_cont_1_below)[np.asarray(pos_1_below)]))
                mask_slit_1_below.append(0)

            xs_cont_slit=np.asarray(xs_cont_1_below)[np.asarray(pos_1_below_1)]
            zs_cont_slit=np.asarray(zs_cont_1_below)[np.asarray(pos_1_below_1)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2)):
                r_slit_1_below.append(-j)
                v_slit_1_below.append(10)
                mask_slit_1_below.append(1)

            else:
                r_slit_1_below.append(-j)
                v_slit_1_below.append(np.mean(np.asarray(vel_cont_1_below)[np.asarray(pos_1_below_1)]))
                mask_slit_1_below.append(0)

        '''Now PA 2 which is greater than 90'''
        width=pixelscale_2*j

        z_2_perp_upper=[-x*(1/np.tan(PA_obs_2))+(width+1.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2]
        z_2_perp_lower=[-x*(1/np.tan(PA_obs_2))+(width+0.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2]
        z_2_perpl_upper=[-x*(1/np.tan(PA_obs_2))-(width+0.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2]
        z_2_perpl_lower=[-x*(1/np.tan(PA_obs_2))-(width+1.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2]
        
        
        if PA2 > 90:
            pos_1_2=np.logical_and((np.asarray(zs_cont_2) < np.asarray(z_2_perp_upper)), (np.asarray(zs_cont_2) > np.asarray(z_2_perp_lower)))
            pos_1_3=np.logical_and((np.asarray(zs_cont_2) < np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2)> np.asarray(z_2_perpl_lower)))

            pos_1_2_center=np.logical_and((np.asarray(zs_cont_2) > np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2) < np.asarray(z_2_perp_lower)))
        else:
            pos_1_2=np.logical_and((np.asarray(zs_cont_2) > np.asarray(z_2_perp_upper)), (np.asarray(zs_cont_2) < np.asarray(z_2_perp_lower)))
            pos_1_3=np.logical_and((np.asarray(zs_cont_2) > np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2)< np.asarray(z_2_perpl_lower)))

            pos_1_2_center=np.logical_and((np.asarray(zs_cont_2) < np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2) > np.asarray(z_2_perp_lower)))
        #
        
        #if j==0:
            #plt.plot(xs_cont_2, z_2_perp_upper, label='z_1_perp_upper',lw=4, color='black')
            #plt.plot(xs_cont_2, z_2_perp_lower, label='z_2_perp_lower',lw=4,ls='--', color='red')
            #plt.plot(xs_cont_2, z_2_perpl_upper, label='z_2_perpl_upper',lw=4, color='pink')
            #plt.show()
        if j==0:
            xs_cont_slit_middle=np.asarray(xs_cont_2)[np.asarray(pos_1_2_center)]
            zs_cont_slit_middle=np.asarray(zs_cont_2)[np.asarray(pos_1_2_center)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit_middle)**2+np.mean(zs_cont_slit_middle)**2))==True:

                r_slit_2_above.append(j)
                v_slit_2_above.append(10)
                mask_slit_2_above.append(1)

            else:
                #so this is the positions to the North so we'll make them negative
                r_slit_2_above.append(j)
                v_slit_2_above.append(np.mean(np.asarray(vel_cont_2)[np.asarray(pos_1_2_center)]))
                mask_slit_2_above.append(0)
        else:

            xs_cont_slit=np.asarray(xs_cont_2)[np.asarray(pos_1_2)]
            zs_cont_slit=np.asarray(zs_cont_2)[np.asarray(pos_1_2)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))==True:
                r_slit_2_above.append(j)
                v_slit_2_above.append(10)
                mask_slit_2_above.append(1)

            else:
                r_slit_2_above.append(j)
                v_slit_2_above.append(np.mean(np.asarray(vel_cont_2)[np.asarray(pos_1_2)]))
                mask_slit_2_above.append(0)

            xs_cont_slit=np.asarray(xs_cont_2)[np.asarray(pos_1_3)]
            zs_cont_slit=np.asarray(zs_cont_2)[np.asarray(pos_1_3)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))==True:

                r_slit_2_above.append(-j)
                v_slit_2_above.append(10)
                mask_slit_2_above.append(1)

            else:
                r_slit_2_above.append(-j)
                v_slit_2_above.append(np.mean(np.asarray(vel_cont_2)[np.asarray(pos_1_3)]))
                mask_slit_2_above.append(0)
            #math.isnan(np.mean(np.asarray(vel_cont_2)[np.asarray(pos_1_3)])):



        z_2_perp_upper=[-x*(1/np.tan(PA_obs_2))+(width+1.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2_below]
        z_2_perp_lower=[-x*(1/np.tan(PA_obs_2))+(width+0.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2_below]
        z_2_perpl_upper=[-x*(1/np.tan(PA_obs_2))-(width+0.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2_below]
        z_2_perpl_lower=[-x*(1/np.tan(PA_obs_2))-(width+1.5*pixelscale_2)/np.cos(PA_obs_2+additive) for x in xs_cont_2_below]

        
        if PA2 >  90:
            pos_1_4=np.logical_and((np.asarray(zs_cont_2_below) < np.asarray(z_2_perp_upper)), (np.asarray(zs_cont_2_below) > np.asarray(z_2_perp_lower)))
            pos_1_5=np.logical_and((np.asarray(zs_cont_2_below) < np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2_below) > np.asarray(z_2_perpl_lower)))
            pos_1_4_center=np.logical_and((np.asarray(zs_cont_2_below) > np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2_below) < np.asarray(z_2_perp_lower)))
        else:
            pos_1_4=np.logical_and((np.asarray(zs_cont_2_below) > np.asarray(z_2_perp_upper)), (np.asarray(zs_cont_2_below) < np.asarray(z_2_perp_lower)))
            pos_1_5=np.logical_and((np.asarray(zs_cont_2_below) > np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2_below) < np.asarray(z_2_perpl_lower)))
            pos_1_4_center=np.logical_and((np.asarray(zs_cont_2_below) < np.asarray(z_2_perpl_upper)), (np.asarray(zs_cont_2_below) > np.asarray(z_2_perp_lower)))


        if j==0:
            xs_cont_slit_middle=np.asarray(xs_cont_2_below)[np.asarray(pos_1_4_center)]
            zs_cont_slit_middle=np.asarray(zs_cont_2_below)[np.asarray(pos_1_4_center)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit_middle)**2+np.mean(zs_cont_slit_middle)**2))==True:

                r_slit_2_below.append(j)
                v_slit_2_below.append(10)
                mask_slit_2_below.append(1)
            else:
                #so this is the positions to the North so we'll make them negative
                r_slit_2_below.append(j)
                v_slit_2_below.append(np.mean(np.asarray(vel_cont_2_below)[np.asarray(pos_1_4_center)]))
                mask_slit_2_below.append(0)
        else:

            xs_cont_slit=np.asarray(xs_cont_2_below)[np.asarray(pos_1_4)]
            zs_cont_slit=np.asarray(zs_cont_2_below)[np.asarray(pos_1_4)]
            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))==True:
                r_slit_2_below.append(j)
                v_slit_2_below.append(10)
                mask_slit_2_below.append(1)

            else:
                r_slit_2_below.append(j)
                v_slit_2_below.append(np.mean(np.asarray(vel_cont_2_below)[np.asarray(pos_1_4)]))
                mask_slit_2_below.append(0)

            xs_cont_slit=np.asarray(xs_cont_2_below)[np.asarray(pos_1_5)]
            zs_cont_slit=np.asarray(zs_cont_2_below)[np.asarray(pos_1_5)]

            if math.isnan(np.sqrt(np.mean(xs_cont_slit)**2+np.mean(zs_cont_slit)**2))==True:

                r_slit_2_below.append(-j)
                v_slit_2_below.append(10)
                mask_slit_2_below.append(1)

            else:
                r_slit_2_below.append(-j)
                v_slit_2_below.append(np.mean(np.asarray(vel_cont_2_below)[np.asarray(pos_1_5)]))
                mask_slit_2_below.append(0)



    inds = np.array(r_slit_1_above).argsort()
    sortedx_1_above = np.array(r_slit_1_above)[inds]#-minus_one
    sortedv_1_above = np.array(v_slit_1_above)[inds]

    sorted_mask = np.array(mask_slit_1_above)[inds]

    m_sortedx_1_above = ma.masked_array(sortedx_1_above, mask=sorted_mask)
    m_sortedv_1_above = ma.masked_array(sortedv_1_above, mask=sorted_mask)

    inds = np.array(r_slit_1_below).argsort()
    sortedx_1_below = np.array(r_slit_1_below)[inds]
    sortedv_1_below = np.array(v_slit_1_below)[inds]

    sorted_mask = np.array(mask_slit_1_below)[inds]

    m_sortedx_1_below = ma.masked_array(sortedx_1_below, mask=sorted_mask)
    m_sortedv_1_below = ma.masked_array(sortedv_1_below, mask=sorted_mask)

    inds = np.array(r_slit_2_above).argsort()
    sortedx_2_above = np.array(r_slit_2_above)[inds]
    sortedv_2_above = np.array(v_slit_2_above)[inds]

    sorted_mask = np.array(mask_slit_2_above)[inds]

    m_sortedx_2_above = ma.masked_array(sortedx_2_above, mask=sorted_mask)
    m_sortedv_2_above = ma.masked_array(sortedv_2_above, mask=sorted_mask)

    inds = np.array(r_slit_2_below).argsort()
    sortedx_2_below = np.array(r_slit_2_below)[inds]
    sortedv_2_below = np.array(v_slit_2_below)[inds]

    sorted_mask = np.array(mask_slit_2_below)[inds]

    m_sortedx_2_below = ma.masked_array(sortedx_2_below, mask=sorted_mask)
    m_sortedv_2_below = ma.masked_array(sortedv_2_below, mask=sorted_mask)



    '''Okay so we need to cut this according to the xs_data_concat array'''
    first=xs_data[lower_lim:upper_lim][0]
    last=xs_data[lower_lim:upper_lim][-1]
    for z in range(len(m_sortedx_1_below)):
        if sortedx_1_below[z]==first:
            i_b_first=z
        if sortedx_1_below[z]==last:
            i_b_last=z+1

    for z in range(len(m_sortedx_1_above)):
        if sortedx_1_above[z]==first:
            i_a_first=z
        if sortedx_1_above[z]==last:
            i_a_last=z+1


    xs_fit_concat=np.ma.concatenate([m_sortedx_1_below[i_b_first:i_b_last],m_sortedx_1_above[i_a_first:i_a_last]])#xs_fit_bin+xs_fit_bin_below
    vs_fit_concat=np.ma.concatenate([m_sortedv_1_below[i_b_first:i_b_last],m_sortedv_1_above[i_a_first:i_a_last]])#vs_fit_bin+vs_fit_bin_below

    first=xs_data_ortho[lower_lim_ortho:upper_lim_ortho][0]
    last=xs_data_ortho[lower_lim_ortho:upper_lim_ortho][-1]

    '''find index where m_sortedx_1_below and same for above have the values of first and last and cut it there'''
    for z in range(len(m_sortedx_2_below)):
        if sortedx_2_below[z]==first:
            i_b_first=z
        if sortedx_2_below[z]==last:
            i_b_last=z+1

    for z in range(len(m_sortedx_2_above)):
        if sortedx_2_above[z]==first:
            i_a_first=z
        if sortedx_2_above[z]==last:
            i_a_last=z+1




    xs_fit_concat_ortho=np.ma.concatenate([m_sortedx_2_below[i_b_first:i_b_last],m_sortedx_2_above[i_a_first:i_a_last]])#xs_fit_ortho_bin+xs_fit_ortho_bin_below
    vs_fit_concat_ortho=np.ma.concatenate([m_sortedv_2_below[i_b_first:i_b_last],m_sortedv_2_above[i_a_first:i_a_last]])#vs_fit_ortho_bin+vs_fit_ortho_bin_below


    xs_data_concat_arc=[x*pixelscale for x in xs_data_concat]
    xs_data_concat_ortho_arc=[x*pixelscale for x in xs_data_concat_ortho]
    xs_fit_concat_arc=[x*pixelscale for x in xs_fit_concat]
    xs_fit_concat_ortho_arc=[x*pixelscale for x in xs_fit_concat_ortho]

 
    
    
    resid_alt_1=[]
    resid_alt_2=[]

    sum_masked=[]
    for i in range(len(vs_fit_concat)):
        if vs_fit_concat[i] is ma.masked:
            sum_masked.append(1)
        else:
            resid_alt_1.append(-0.5*(((vs_fit_concat[i]-vel_data_concat[i])**2)/(np.sqrt((np.sqrt(error_data_concat[0][i]**2+error_data_concat[1][i]**2))**2+2.5**2))**2))
    for i in range(len(vs_fit_concat_ortho)):
        if vs_fit_concat_ortho[i] is ma.masked:
            sum_masked.append(1)
        else:
            resid_alt_2.append(-0.5*(((vs_fit_concat_ortho[i]-vel_data_concat_ortho[i])**2)/(np.sqrt((np.sqrt(error_data_concat_ortho[0][i]**2+error_data_concat_ortho[1][i]**2))**2+2.5**2))**2))

    masked_elements_sum=2*np.sum(sum_masked)

    resid_alt_concat=resid_alt_1+resid_alt_2
    #resid_concat=resid_1+resid_2




    n=2*(upper_lim-lower_lim+upper_lim_ortho-lower_lim_ortho)
    k=5



    return np.sum(resid_alt_concat)/(n-k-1-masked_elements_sum)









# Now, this part calls whatever galaxy you are running
names=['Akira']
for z in range(len(names)):
    name=names[z]
    
    print(name)
    
    minus_one=0
    
    if name=='Akira':
        PA_1 = 'akira_slit1_input.txt'
        PA_2 = 'akira_slit2_input.txt'
        PA1 = 40
        PA2 = 110
        lower_lim= 0
        upper_lim = 40
        lower_lim_ortho = 0
        upper_lim_ortho = 32
        pixelscale = 0.15
        pixelscale_1 = pixelscale
        pixelscale_2 = pixelscale
        slitwidth = 0.75
        result_s = 45, 55, 10, 30, 310
        
        minus_one=1
        minus_pa_1 = 2
    
    

    with open(PA_1, 'r') as f:
        data = f.readlines()
        title=[]
        col_1=[]
        col_2=[]
        col_3=[]
        col_1_1=[]
        col_1_2=[]
        col_1_3=[]
        col_2_1=[]
        col_2_2=[]
        col_2_3=[]

        for line in data:
            words = line.split()
            if words[0]=='Vel':
                title.append(words[0])
                col_1.append(words[1])
                col_2.append(words[3])
                col_3.append(words[4])
            if words[0]=='Veltwoone':
                col_1_1.append(words[1])
                col_1_2.append(words[3])
                col_1_3.append(words[4])
            if words[0]=='Veltwotwo':
                col_2_1.append(words[1])
                col_2_2.append(words[3])
                col_2_3.append(words[4])

    vel_1_gauss=[]
    vel_1_error_up=[]
    vel_1_error_down=[]
    vel_1_error=[]

    vel_2_1_gauss=[]
    vel_2_1_error_up=[]
    vel_2_1_error_down=[]
    vel_2_1_error=[]

    vel_2_2_gauss=[]
    vel_2_2_error_up=[]
    vel_2_2_error_down=[]
    vel_2_2_error=[]

    for j in range(len(col_1_1)):
        if j < 40:
            vel_1_gauss.append(col_1[j])
            vel_1_error_up.append(float(col_2[j])-float(col_1[j]))
            vel_1_error_down.append(float(col_1[j])-float(col_3[j]))
            vel_1_error.append([float(col_2[j])-float(col_1[j]),float(col_1[j])-float(col_3[j])])

            vel_2_1_gauss.append(col_1_1[j])
            if name=='Akira':
                vel_2_1_error_up.append(float(col_1_2[j]))
                vel_2_1_error_down.append(float(col_1_1[j]))
                vel_2_1_error.append([float(col_1_2[j]),float(col_1_3[j])])
                vel_2_2_gauss.append(col_2_1[j])
                vel_2_2_error_up.append(float(col_2_2[j]))
                vel_2_2_error_down.append(float(col_2_1[j]))
                vel_2_2_error.append([float(col_2_2[j]),float(col_2_3[j])])       
            else:

                vel_2_1_error_up.append(float(col_1_2[j])-float(col_1_1[j]))
                vel_2_1_error_down.append(float(col_1_1[j])-float(col_1_3[j]))
                vel_2_1_error.append([float(col_1_2[j])-float(col_1_1[j]),float(col_1_1[j])-float(col_1_3[j])])

                vel_2_2_gauss.append(col_2_1[j])
                vel_2_2_error_up.append(float(col_2_2[j])-float(col_2_1[j]))
                vel_2_2_error_down.append(float(col_2_1[j])-float(col_2_3[j]))
                vel_2_2_error.append([float(col_2_2[j])-float(col_2_1[j]),float(col_2_1[j])-float(col_2_3[j])])                                                                                                                                                                                                


    
    with open(PA_2, 'r') as f:
        data = f.readlines()
        title=[]
        col_1=[]
        col_2=[]
        col_3=[]
        col_1_1=[]
        col_1_2=[]
        col_1_3=[]
        col_2_1=[]
        col_2_2=[]
        col_2_3=[]

        for line in data:
            words = line.split()
            if words[0]=='Vel':
                title.append(words[0])
                col_1.append(words[1])
                col_2.append(words[3])
                col_3.append(words[4])
            if words[0]=='Veltwoone':
                col_1_1.append(words[1])
                col_1_2.append(words[3])
                col_1_3.append(words[4])
            if words[0]=='Veltwotwo':
                col_2_1.append(words[1])
                col_2_2.append(words[3])
                col_2_3.append(words[4])


    vel_1_gauss_ortho=[]
    vel_1_error_ortho_up=[]
    vel_1_error_ortho_down=[]
    vel_1_error_ortho=[]

    vel_2_1_gauss_ortho=[]
    vel_2_1_error_ortho_up=[]
    vel_2_1_error_ortho_down=[]
    vel_2_1_error_ortho=[]

    vel_2_2_gauss_ortho=[]
    vel_2_2_error_ortho_up=[]
    vel_2_2_error_ortho_down=[]
    vel_2_2_error_ortho=[]

    for j in range(len(col_1_1)):
        if j < 40:
            vel_1_gauss_ortho.append(col_1[j])
            vel_1_error_ortho_up.append(float(col_2[j])-float(col_1[j]))
            vel_1_error_ortho_down.append(float(col_1[j])-float(col_3[j]))
            vel_1_error_ortho.append([float(col_2[j])-float(col_1[j]),float(col_1[j])-float(col_3[j])])


            if name=='Akira':
                vel_2_1_gauss_ortho.append(col_1_1[j])
                vel_2_1_error_ortho_up.append(float(col_1_2[j]))
                vel_2_1_error_ortho_down.append(float(col_1_3[j]))
                vel_2_1_error_ortho.append([float(col_1_2[j]),float(col_1_3[j])])


                vel_2_2_gauss_ortho.append(col_2_1[j])
                vel_2_2_error_ortho_up.append(float(col_2_2[j]))
                vel_2_2_error_ortho_down.append(float(col_2_3[j]))
                vel_2_2_error_ortho.append([float(col_2_2[j]),float(col_2_3[j])])

            else:
                vel_2_1_gauss_ortho.append(col_1_1[j])
                vel_2_1_error_ortho_up.append(float(col_1_2[j])-float(col_1_1[j]))
                vel_2_1_error_ortho_down.append(float(col_1_1[j])-float(col_1_3[j]))
                vel_2_1_error_ortho.append([float(col_1_2[j])-float(col_1_1[j]),float(col_1_1[j])-float(col_1_3[j])])

                vel_2_2_gauss_ortho.append(col_2_1[j])
                vel_2_2_error_ortho_up.append(float(col_2_2[j])-float(col_2_1[j]))
                vel_2_2_error_ortho_down.append(float(col_2_1[j])-float(col_2_3[j]))
                vel_2_2_error_ortho.append([float(col_2_2[j])-float(col_2_1[j]),float(col_2_1[j])-float(col_2_3[j])])
            





    xs_data=np.linspace(0, len(vel_1_gauss)-1, len(vel_1_gauss))

    xs_data=xs_data-len(xs_data)/2 - minus_pa_1


    xs_data=[float(y) for y in xs_data]
    vel_1=[float(y) for y in vel_1_gauss]
    vel_2_1=[float(y) for y in vel_2_1_gauss]
    vel_2_2=[float(y) for y in vel_2_2_gauss]
    #vel_1_error=[float(y) for y in vel_1_error]

    xs_data=[float(y) for y in xs_data]
    vel_1_ortho=[float(y) for y in vel_1_gauss_ortho]
    vel_2_1_ortho=[float(y) for y in vel_2_1_gauss_ortho]
    vel_2_2_ortho=[float(y) for y in vel_2_2_gauss_ortho]

    #vel_1_error_ortho=[float(y) for y in vel_1_error_ortho]


    # So at this point, ALL rows of the data table are imported, after this then they are restricted:

    if np.mean(vel_2_1[lower_lim:upper_lim])< np.mean(vel_2_1[lower_lim:upper_lim]):
        #this means that the first set is the above side
        v_1_above=vel_2_1[lower_lim:upper_lim]
        v_1_below=vel_2_2[lower_lim:upper_lim]
        vel_data_concat=vel_2_1[lower_lim:upper_lim]+vel_2_2[lower_lim:upper_lim]
        vel_error_up_concat=vel_2_1_error_up[lower_lim:upper_lim]+vel_2_2_error_up[lower_lim:upper_lim]
        vel_error_down_concat=vel_2_1_error_down[lower_lim:upper_lim]+vel_2_2_error_down[lower_lim:upper_lim]
        #error_data_concat=vel_2_1_error[lower_lim:upper_lim]+vel_2_2_error[lower_lim:upper_lim]
    else:
        v_1_above=vel_2_2[lower_lim:upper_lim]
        v_1_below=vel_2_1[lower_lim:upper_lim]
        vel_data_concat=vel_2_2[lower_lim:upper_lim]+vel_2_1[lower_lim:upper_lim]
        vel_error_up_concat=vel_2_2_error_up[lower_lim:upper_lim]+vel_2_1_error_up[lower_lim:upper_lim]
        vel_error_down_concat=vel_2_2_error_down[lower_lim:upper_lim]+vel_2_1_error_down[lower_lim:upper_lim]
        #error_data_concat=vel_2_2_error[lower_lim:upper_lim]+vel_2_1_error[lower_lim:upper_lim]

    


    if np.mean(vel_2_1_ortho[lower_lim_ortho:upper_lim_ortho])< np.mean(vel_2_1_ortho[lower_lim_ortho:upper_lim_ortho]):
        #this means that the first set is the above side
        v_2_above=vel_2_1_ortho[lower_lim_ortho:upper_lim_ortho]
        v_2_below=vel_2_2_ortho[lower_lim_ortho:upper_lim_ortho]
        vel_data_concat_ortho=vel_2_1_ortho[lower_lim_ortho:upper_lim_ortho]+vel_2_2_ortho[lower_lim_ortho:upper_lim_ortho]
        vel_error_ortho_up_concat=vel_2_1_error_ortho_up[lower_lim_ortho:upper_lim_ortho]+vel_2_2_error_ortho_up[lower_lim_ortho:upper_lim_ortho]
        vel_error_ortho_down_concat=vel_2_1_error_ortho_down[lower_lim_ortho:upper_lim_ortho]+vel_2_2_error_ortho_down[lower_lim_ortho:upper_lim_ortho]
        #error_data_concat_ortho=vel_2_1_ortho_error[lower_lim_ortho:upper_lim_ortho]+vel_2_2_ortho_error[lower_lim_ortho:upper_lim_ortho]
    else:
        v_2_above=vel_2_2_ortho[lower_lim_ortho:upper_lim_ortho]
        v_2_below=vel_2_1_ortho[lower_lim_ortho:upper_lim_ortho]
        vel_data_concat_ortho=vel_2_2_ortho[lower_lim_ortho:upper_lim_ortho]+vel_2_1_ortho[lower_lim_ortho:upper_lim_ortho]
        vel_error_ortho_up_concat=vel_2_2_error_ortho_up[lower_lim_ortho:upper_lim_ortho]+vel_2_1_error_ortho_up[lower_lim_ortho:upper_lim_ortho]
        vel_error_ortho_down_concat=vel_2_2_error_ortho_down[lower_lim_ortho:upper_lim_ortho]+vel_2_1_error_ortho_down[lower_lim_ortho:upper_lim_ortho]

    xs_data_cut=xs_data[lower_lim:upper_lim]
    xs_data_concat=xs_data_cut+xs_data_cut
    xs_data_concat=[int(x) for x in xs_data_concat]



    xs_data_ortho=np.linspace(0, len(vel_1_gauss_ortho)-1, len(vel_1_gauss_ortho))

    xs_data_ortho=xs_data_ortho-len(xs_data_ortho)/2+0.5


    xs_data_ortho=[float(y) for y in xs_data_ortho]
    print('xs_data_ortho', xs_data_ortho)

    xs_data_cut_ortho= xs_data_ortho[lower_lim_ortho:upper_lim_ortho]
    xs_data_concat_ortho=xs_data_cut_ortho+xs_data_cut_ortho
    xs_data_concat_ortho=[int(x) for x in xs_data_concat_ortho]

    print('cut', xs_data_cut_ortho, lower_lim_ortho, upper_lim_ortho)
    print(xs_data_concat_ortho)



    xs_data_concat_arc=[pixelscale_1*y for y in xs_data_concat]
    xs_data_concat_ortho_arc=[pixelscale_2*y for y in xs_data_concat_ortho]



    assy_error=[vel_error_down_concat,vel_error_up_concat]
    assy_error_ortho=[vel_error_ortho_down_concat,vel_error_ortho_up_concat]

    error_data_concat=assy_error
    error_data_concat_ortho=assy_error_ortho


    start=datetime.datetime.now()
    z_prelim=10,30,2,40,500

    print('xs_data', xs_data_concat)
    print('vel_data', vel_data_concat)
    print('xs_data_ortho', xs_data_concat_ortho)
    print('vel_data_ortho', vel_data_concat_ortho)
    print(len(xs_data_concat), len(vel_data_concat),len(xs_data_concat_ortho), len(vel_data_concat_ortho))
    
    import seaborn as sns
    sns.set_style("whitegrid")
    plt.clf()
    fig = plt.figure()
    ax0 = fig.add_subplot(211)
    ax0.scatter(xs_data_concat, vel_data_concat, label='PA '+str(PA1), color='#905569')
    ax0.errorbar(xs_data_concat, vel_data_concat, yerr=vel_error_down_concat, capsize=5, 
                 color='#905569',ls='None')
    ax0.legend()
    ax0.set_ylim([-500,500])

    ax1 = fig.add_subplot(212)
    ax1.scatter(xs_data_concat_ortho, vel_data_concat_ortho, label='PA '+str(PA2), color='#E09F6A')
    ax1.errorbar(xs_data_concat_ortho, vel_data_concat_ortho, yerr=vel_error_ortho_down_concat, capsize=5, 
                 ls='None', color='#E09F6A')
    ax1.legend()
    ax1.set_ylim([-500,500])
    plt.savefig('data.png')



    # Maybe also plot the slits against the 2D map:
    xs_2D_concat_1 = [pixelscale_1*np.cos(math.radians(PA1))*x+(pixelscale_1*np.cos(math.radians(PA1)))/2 for x in xs_data_concat[:-int(len(xs_data_concat)/2)]]
    # pixel scale is "/pix, so this whole plot will be in arcsec
    ys_2D_concat_1 = [pixelscale_1*np.sin(math.radians(PA1))*x-(pixelscale_1*np.cos(math.radians(PA1)))/2 for x in xs_data_concat[:-int(len(xs_data_concat)/2)]]
    xs_2D_concat_2 = [pixelscale_1*np.cos(math.radians(PA1))*x-(pixelscale_1*np.cos(math.radians(PA1)))/2 for x in xs_data_concat[-int(len(xs_data_concat)/2):]]
    ys_2D_concat_2 = [pixelscale_1*np.sin(math.radians(PA1))*x+(pixelscale_1*np.cos(math.radians(PA1)))/2 for x in xs_data_concat[-int(len(xs_data_concat)/2):]]
    


    xs_2D_concat_ortho_1 = [pixelscale_2*np.cos(math.radians(PA2))*x+(pixelscale_2*np.sin(math.radians(PA2)))/2 for x in xs_data_concat_ortho[:-int(len(xs_data_concat_ortho)/2)]]
    # pixel scale is "/pix, so this whole plot will be in arcsec
    ys_2D_concat_ortho_1 = [pixelscale_2*np.sin(math.radians(PA2))*x for x in xs_data_concat_ortho[:-int(len(xs_data_concat_ortho)/2)]]
    xs_2D_concat_ortho_2 = [pixelscale_2*np.cos(math.radians(PA2))*x-(pixelscale_2*np.sin(math.radians(PA2)))/2 for x in xs_data_concat_ortho[-int(len(xs_data_concat_ortho)/2):]]
    ys_2D_concat_ortho_2 = [pixelscale_2*np.sin(math.radians(PA2))*x for x in xs_data_concat_ortho[-int(len(xs_data_concat_ortho)/2):]]

    plt.clf()
    plt.scatter(xs_2D_concat_1, ys_2D_concat_1, c=vel_data_concat[:-int(len(xs_data_concat)/2)], cmap='RdBu_r', vmax=300, vmin=-300)
    plt.scatter(xs_2D_concat_2, ys_2D_concat_2, c=vel_data_concat[-int(len(xs_data_concat)/2):], cmap='RdBu_r', vmax=300, vmin=-300)
    plt.scatter(xs_2D_concat_ortho_1, ys_2D_concat_ortho_1, c=vel_data_concat_ortho[:-int(len(xs_data_concat_ortho)/2)], cmap='RdBu_r', vmax=300, vmin=-300)
    plt.scatter(xs_2D_concat_ortho_2, ys_2D_concat_ortho_2, c=vel_data_concat_ortho[-int(len(xs_data_concat_ortho)/2):], cmap='RdBu_r', vmax=300, vmin=-300)
    
    plt.colorbar()
    plt.savefig('physical.png')
    


    # At this point, we have prepared the data
    # Try just plotting this against the model
    
    input = result_s
    xs_model, ys_model, xs_model_ortho, ys_model_ortho = diagnostic_bicone_plot(input, xs_data_concat, xs_data_concat_ortho)
    print(np.shape(xs_model))
    plt.clf()
    fig = plt.figure()
    ax0 = fig.add_subplot(211)
    
    ax0.scatter(xs_model, ys_model, label='Model PA1', color='orange')
    print('value', ys_model)
    
    print('error down',vel_error_down_concat)
    print('error up', vel_error_up_concat)
    #plt.fill_between(xs_model, ys_model - vel_error_down_concat, ys_model + vel_error_up_concat, alpha=0.5, color='orange') 
    #ax0.errorbar(xs_model, ys_model, yerr=[vel_error_down_concat, vel_error_up_concat], fmt='ko', label='Model PA1',capsize=5)
    
    
    ax1 = fig.add_subplot(212)
    #ax1.scatter(xs_model_ortho, ys_model_ortho, label='Model PA2', color='red')
    ax0.scatter(xs_data_concat, vel_data_concat, label='PA1', color='green')
    ax0.legend(loc='upper right')
    #plt.errorbars(xs_data_concat, vel_data_concat, ve
    ax1.scatter(xs_data_concat_ortho, vel_data_concat_ortho, label='PA2',color='green')
    print(len(xs_model_ortho), len(ys_model_ortho), len(vel_error_ortho_down_concat))
    #ax1.errorbar(xs_model_ortho, ys_model_ortho, 
    #             yerr=[vel_error_ortho_down_concat, vel_error_ortho_up_concat],fmt='ko', label='Model PA2',
    #             capsize=10)
    ax1.scatter(xs_model_ortho, ys_model_ortho, label='Model PA2', color='orange')
    ax0.annotate(str(input),xy=(0.01,0.9), xycoords='axes fraction') 
    ax0.set_title('PA 1 = '+str(PA1))
    ax0.set_ylim([-500,500])
    ax1.set_title('PA 2 = '+str(PA2))
    #plt.ylim([-500,500])
    #plt.legend()
    plt.legend(loc='upper right')
    plt.ylabel('Velocity [km s$^{-1}$]')
    plt.xlabel('Spatial Position [Rows]')
    plt.tight_layout()
    plt.savefig('testing_Akira.png')
    
    
    #def lnlike_bicone_three_slits(z, PA1, PA2, PA3, pixelscale_1, slitwidth):

    try:
    	lnlike_bicone_three_slits(result_s, PA1, PA2, PA3, pixelscale_1, slitwidth)
    except NameError:
    	lnlike_bicone_three_slits(result_s, PA1, PA2, 888, pixelscale_1, slitwidth)
     
     
    STOP
    lnlike_bicone(result_s)
    
    end=datetime.datetime.now()
    print(end-start)
    print('youve done one run')
    #`~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
    

    def lnprior(z):
        phi_array, theta_array, r_t, half_angle, vel_max_array=z
        if 0 < phi_array < 90 and 0.0 < theta_array < 359.9 and 0 < r_t < 10 and 0 < half_angle < 90 and 10 < vel_max_array < 1000:
            return 0.0
        return -np.inf
    #you also need the log likelihood function:                                                                                                                                                             
    #we already have a log likelihood function                                                                                                                                                              
    def lnprob(z):
        lp=lnprior(z)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike_bicone(z)
    ndim, nwalkers=5,10#was 100
    result=45,18.6,2,63.36,420
    mult=10,20,1,10,100

    #pos is the starting position for the walkers                                                                                                                                                           
    pos=[result+mult*np.random.randn(ndim) for i in range(nwalkers)]
    #print('pos', pos)                                                                                                                                                                                      

    import emcee
    start=datetime.datetime.now()
    #count=0                                                                                                                                                                                                
    sampler=emcee.EnsembleSampler(nwalkers, ndim, lnprob,a=2, threads=48)
    #sampler=emcee.EnsembleSampler(nwalkers, ndim, lnprob,a=1.75, threads=48)
    print('You finished the prep for burn in')
    
    len_chain = 20
    pos, prob, state = sampler.run_mcmc(pos,len_chain, progress=True)#was 200
    end=datetime.datetime.now()
    samples_pre=sampler.chain[:,:,:]
    f_acc_pre=sampler.acceptance_fraction
    #autocorr_pre=sampler.acor
    print('f_acc burn in', f_acc_pre)
    #print('autocorr for burn in', autocorr_pre)
    print('pos after burn in', pos)
    g=open(str(name)+"/bicone/emcee_burn_0_c.txt","w")
    h=open(str(name)+"/bicone/emcee_burn_1_c.txt","w")
    k=open(str(name)+"/bicone/emcee_burn_2_c.txt","w")
    l=open(str(name)+"/bicone/emcee_burn_3_c.txt","w")
    m=open(str(name)+"/bicone/emcee_burn_4_c.txt","w")
    n=open(str(name)+"/bicone/emcee_burn_f_acc_c.txt","w")
    for i in range(nwalkers):
        for j in range(len_chain):
            g.write(str(samples_pre[i][j][0])+'\n')
            h.write(str(samples_pre[i][j][1])+'\n')
            k.write(str(samples_pre[i][j][2])+'\n')
            l.write(str(samples_pre[i][j][3])+'\n')
            m.write(str(samples_pre[i][j][4])+'\n')

    #       #n.write(str(f_acc_pre[i][j])+'\n')                                                                                                                                                             
    g.close()
    h.close()
    k.close()
    l.close()
    m.close()
    
    
    plt.clf()
    fig, axes = plt.subplots(5, figsize=(15, 7), sharex=True)
    samples = sampler.get_chain()
    labels = ["inc", "PA", "r","theta", "v"]
    for i in range(ndim):
    	ax = axes[i]
    	ax.plot(samples[:, :, i], "k", alpha=0.3)
    	ax.set_xlim(0, len(samples))
    	ax.set_ylabel(labels[i])
    	ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")
    plt.savefig('walkers_bicone.png')


    print(np.shape(samples_pre))
    mean_0 = np.mean(samples_pre[:][len_chain][0])
    mean_1 = np.mean(samples_pre[:][len_chain][1])
    mean_2 = np.mean(samples_pre[:][len_chain][2])
    mean_3 = np.mean(samples_pre[:][len_chain][3])
    mean_4 = np.mean(samples_pre[:][len_chain][4])
    print('mean variables', mean_0, mean_1, mean_2, mean_3, mean_4)
    

    STOP

    for i in range(len(f_acc_pre)):
        n.write(str(f_acc_pre[i])+'\n')
    n.close()
    sampler.reset()

    def lnprior(z):
        phi_array, theta_array, r_t, half_angle,half_angle_2, vel_max_array=z
        if 0 < phi_array < 90 and 0.0 < theta_array < 359.9 and 0 < r_t < 10 and 0 < half_angle < 90 and half_angle < half_angle_2 < 90 and 100 < vel_max_array < 1000:
            return 0.0
        return -np.inf
    #you also need the log likelihood function:                                                                                                                                                             
    #we already have a log likelihood function                                                                                                                                                              
    def lnprob(z):
        lp=lnprior(z)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike_cocone(z)
    ndim, nwalkers=6,100
    result=45,18.6,2,63.36,72,420
    mult=10,20,1,10,10,100

    #pos is the starting position for the walkers                                                                                                                                                           
    pos=[result+mult*np.random.randn(ndim) for i in range(nwalkers)]
    #print('pos', pos)                                                                                                                                                                                      

    import emcee
    start=datetime.datetime.now()
    #count=0                                                                                                                                                                                                
    sampler=emcee.EnsembleSampler(nwalkers, ndim, lnprob,a=2, threads=48)
    #sampler=emcee.EnsembleSampler(nwalkers, ndim, lnprob,a=1.75, threads=48)
    print('You finished the prep for burn in')
    pos, prob, state = sampler.run_mcmc(pos,200)#, storechain=True)
    end=datetime.datetime.now()
    samples_pre=sampler.chain[:,:,:]
    f_acc_pre=sampler.acceptance_fraction
    #autocorr_pre=sampler.acor
    print('f_acc burn in', f_acc_pre)
    #print('autocorr for burn in', autocorr_pre)
    print('pos after burn in', pos)
    g=open(str(name)+"/cocone/emcee_burn_0_c.txt","w")
    h=open(str(name)+"/cocone/emcee_burn_1_c.txt","w")
    k=open(str(name)+"/cocone/emcee_burn_2_c.txt","w")
    l=open(str(name)+"/cocone/emcee_burn_3_c.txt","w")
    m=open(str(name)+"/cocone/emcee_burn_4_c.txt","w")
    o=open(str(name)+"/cocone/emcee_burn_5_c.txt","w")
    n=open(str(name)+"/cocone/emcee_burn_f_acc_c.txt","w")
    for i in range(nwalkers):
        for j in range(200):
            g.write(str(samples_pre[i][j][0])+'\n')
            h.write(str(samples_pre[i][j][1])+'\n')
            k.write(str(samples_pre[i][j][2])+'\n')
            l.write(str(samples_pre[i][j][3])+'\n')
            m.write(str(samples_pre[i][j][4])+'\n')
            o.write(str(samples_pre[i][j][5])+'\n')
    #       #n.write(str(f_acc_pre[i][j])+'\n')                                                                                                                                                             
    g.close()
    h.close()
    k.close()
    l.close()
    m.close()
    o.close()

    for i in range(len(f_acc_pre)):
        n.write(str(f_acc_pre[i])+'\n')
    n.close()
    sampler.reset()


    def lnprior(z):
        phi_array, theta_array, r_t, half_angle,half_angle_2, vel_max_array=z
        if 0 < phi_array < 90 and 0.0 < theta_array < 359.9 and 0 < r_t < 10 and 0 < half_angle < 90 and half_angle < half_angle_2 < 90 and 100 < vel_max_array < 1000:
            return 0.0
        return -np.inf
    #you also need the log likelihood function:                                                                                                                                                             
    #we already have a log likelihood function                                                                                                                                                              
    def lnprob(z):
        lp=lnprior(z)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike_nesting(z)
    ndim, nwalkers=6,10# was 100
    result=45,18.6,2,63.36,72,420
    mult=10,20,1,10,10,100

    #pos is the starting position for the walkers                                                                                                                                                           
    pos=[result+mult*np.random.randn(ndim) for i in range(nwalkers)]
    #print('pos', pos)                                                                                                                                                                                      

    import emcee
    start=datetime.datetime.now()
    #count=0                                                                                                                                                                                                
    sampler=emcee.EnsembleSampler(nwalkers, ndim, lnprob,a=2, threads=48)
    #sampler=emcee.EnsembleSampler(nwalkers, ndim, lnprob,a=1.75, threads=48)
    print('You finished the prep for burn in')
    pos, prob, state = sampler.run_mcmc(pos,20,progress=True)#was 200 steps
    end=datetime.datetime.now()
    samples_pre=sampler.chain[:,:,:]
    f_acc_pre=sampler.acceptance_fraction
    #autocorr_pre=sampler.acor
    print('f_acc burn in', f_acc_pre)
    #print('autocorr for burn in', autocorr_pre)
    print('pos after burn in', pos)
    g=open(str(name)+"/nesting/emcee_burn_0_c.txt","w")
    h=open(str(name)+"/nesting/emcee_burn_1_c.txt","w")
    k=open(str(name)+"/nesting/emcee_burn_2_c.txt","w")
    l=open(str(name)+"/nesting/emcee_burn_3_c.txt","w")
    m=open(str(name)+"/nesting/emcee_burn_4_c.txt","w")
    o=open(str(name)+"/nesting/emcee_burn_5_c.txt","w")
    n=open(str(name)+"/nesting/emcee_burn_f_acc_c.txt","w")
    for i in range(nwalkers):
        for j in range(200):
            g.write(str(samples_pre[i][j][0])+'\n')
            h.write(str(samples_pre[i][j][1])+'\n')
            k.write(str(samples_pre[i][j][2])+'\n')
            l.write(str(samples_pre[i][j][3])+'\n')
            m.write(str(samples_pre[i][j][4])+'\n')
            o.write(str(samples_pre[i][j][5])+'\n')
    #       #n.write(str(f_acc_pre[i][j])+'\n')                                                                                                                                                             
    g.close()
    h.close()
    k.close()
    l.close()
    m.close()
    o.close()

    for i in range(len(f_acc_pre)):
        n.write(str(f_acc_pre[i])+'\n')
    n.close()
    sampler.reset()

    print('done with burnin')


