from __future__ import print_function
import sys
import numpy as np

def read_from_text(file):
	for line in file:
		yield line.strip('\r\n').split('\t')

#shanghai
if 0:
    lat_min=31.05
    lat_max=31.35
    lat_step=0.01
    lat_split=int(np.ceil((lat_max-lat_min)/lat_step))
    lat_max=lat_min+lat_step*lat_split
    
    ratio=np.cos((lat_min+lat_max)*np.pi/360)
    
    lon_min=121.3
    lon_max=121.6
    lon_step=lat_step/ratio
    lon_split=int(np.ceil((lon_max-lon_min)/lon_step))
    lon_max=lon_min+lon_step*lon_split
    
    time_min=1514736000
    time_max=1546272000
    
    time_step=60*60 #time interval
    time_split=(time_max-time_min)/time_step
    #T=24 #per day
    threshold=4*60*60 #hours cut
    #close_time=5*50 #min
    
#shenzhen
if 1:
    lat_min=22.45
    lat_max=22.85
    lat_step=0.01
    lat_split=int(np.ceil((lat_max-lat_min)/lat_step))
    lat_max=lat_min+lat_step*lat_split
    
    ratio=np.cos((lat_min+lat_max)*np.pi/360)
    
    lon_min=113.7
    lon_max=114.4
    lon_step=lat_step/ratio
    lon_split=int(np.ceil((lon_max-lon_min)/lon_step))
    lon_max=lon_min+lon_step*lon_split
    
    time_min=1514736000
    time_max=1546272000
    
    time_step=60*60 #time interval
    time_split=(time_max-time_min)/time_step
    #T=24 #per day
    threshold=4*60*60 #hours cut
    #close_time=5*50 #min

#beijing
if 0:
    lat_min=39.70
    lat_max=40.20
    lat_step=0.01
    lat_split=int(np.ceil((lat_max-lat_min)/lat_step))
    lat_max=lat_min+lat_step*lat_split
    
    ratio=np.cos((lat_min+lat_max)*np.pi/360)
    
    lon_min=116.10
    lon_max=116.70
    lon_step=lat_step/ratio
    lon_split=int(np.ceil((lon_max-lon_min)/lon_step))
    lon_max=lon_min+lon_step*lon_split
    
    time_min=1514736000
    time_max=1546272000
    
    time_step=60*60 #time interval
    time_split=(time_max-time_min)/time_step
    #T=24 #per day
    threshold=4*60*60 #hours cut
    #close_time=5*50 #min



N=lon_split*lat_split
M=np.arange(N,dtype=np.int32).reshape([lon_split,lat_split])

def get_traj_from_line(line):
    points = line[1].split(";")
    traj_len = len(points)
    time_traj = np.zeros(traj_len,dtype=np.int64)
    lon_traj = np.zeros(traj_len,dtype=np.float64)
    lat_traj = np.zeros(traj_len,dtype=np.float64)
    
    point_index = 0
    for p in points:
        lon, lat, time = p.split(",")
        lon, lat, time = float(lon), float(lat), int(time)
        time_traj[point_index] = time-time_min
        lon_traj[point_index] = lon-lon_min
        lat_traj[point_index] = lat-lat_min
        point_index += 1
    
    return time_traj,lon_traj,lat_traj

def get_cutted_trajectory(time_traj,lon_traj,lat_traj):#hours
    if time_traj.shape[0]!=lon_traj.shape[0]:
        print('wrong: lenth unmatch')
        exit(0)
    if time_traj.shape[0]!=lat_traj.shape[0]:
        print('wrong: lenth unmatch')
        exit(0)
        
    cutted_time = []
    cutted_lon = []
    cutted_lat = []
    
    if time_traj.shape[0]>1:
        #cut based on lon lat
        #cut based on time interval
        time_interval = time_traj[1:] - time_traj[:-1]
        time_cut_index = np.where( time_interval > threshold )[0]
        
        cut_index = np.concatenate([np.ones(1,dtype=np.int64)*0,(time_cut_index+np.ones(1,dtype=np.int64)),time_traj.shape],axis=0)
        for i in range(cut_index.shape[0]-1):
            if ( cut_index[i+1] - cut_index[i] ) > 1:
                cutted_time.append(time_traj[cut_index[i]:cut_index[i+1]])
                cutted_lon.append(lon_traj[cut_index[i]:cut_index[i+1]])
                cutted_lat.append(lat_traj[cut_index[i]:cut_index[i+1]])
                
    return cutted_time,cutted_lon,cutted_lat

def get_refined_trajectory(cutted_time_traj,cutted_lon_traj,cutted_lat_traj):
    if cutted_time_traj.shape[0]!=cutted_lon_traj.shape[0]:
        print('wrong: lenth unmatch')
        exit(0)
    if cutted_time_traj.shape[0]!=cutted_lat_traj.shape[0]:
        print('wrong: lenth unmatch')
        exit(0)

    time_min=np.min(cutted_time_traj)
    time_max=np.max(cutted_time_traj)
    time_start=time_min-time_min%time_step+time_step
    time_end=time_max-time_max%time_step
    
    interpolation_num=int((time_end-time_start)/time_step+1)
    
    refined_time_traj=np.arange(time_start,time_end+1,time_step,dtype=np.int64)
    refined_lon_traj=np.zeros(interpolation_num,dtype=np.float64)
    refined_lat_traj=np.zeros(interpolation_num,dtype=np.float64)
    if interpolation_num!=refined_time_traj.shape[0]:
        print('interpolation_num wrong')
        exit(0)
    for i in range(interpolation_num):
        if np.sum(cutted_time_traj==refined_time_traj[i])>0:
            index =np.where(cutted_time_traj==refined_time_traj[i])
            refined_lon_traj[i] =np.mean(cutted_lon_traj[index])
            refined_lat_traj[i] =np.mean(cutted_lat_traj[index])
        else:
            time_less =np.where(cutted_time_traj<refined_time_traj[i])
            time_before =np.max(cutted_time_traj[time_less])
            before_index =np.where(cutted_time_traj==time_before)
            lon_before =np.mean(cutted_lon_traj[before_index])
            lat_before =np.mean(cutted_lat_traj[before_index])
    
            time_larger =np.where(cutted_time_traj>refined_time_traj[i])
            time_later =np.min(cutted_time_traj[time_larger])
            later_index =np.where(cutted_time_traj==time_later)
            lon_later =np.mean(cutted_lon_traj[later_index])
            lat_later =np.mean(cutted_lat_traj[later_index])    
            
            refined_lon_traj[i] = ( lon_before*(time_later-refined_time_traj[i]) + lon_later*(refined_time_traj[i]-time_before) )/(time_later-time_before)
            refined_lat_traj[i] = ( lat_before*(time_later-refined_time_traj[i]) + lat_later*(refined_time_traj[i]-time_before) )/(time_later-time_before)
                
    refined_time_traj=(refined_time_traj/time_step).astype(np.int32)
    refined_lon_traj=np.floor(refined_lon_traj/lon_step).astype(np.int32)
    refined_lat_traj=np.floor(refined_lat_traj/lat_step).astype(np.int32)
    return refined_time_traj,refined_lon_traj,refined_lat_traj

def map_pop(refined_time_traj,refined_lon_traj,refined_lat_traj):
    L=refined_time_traj.shape[0]
    for i in range(L):
        time=refined_time_traj[i]
        lon=refined_lon_traj[i]
        lat=refined_lat_traj[i]
        if (0<=time and time<time_split):
            if (0<=lon and lon<lon_split) and (0<=lat and lat<lat_split):
                location=M[lon,lat]
            else:
                location=N
            print(str(time)+'S'+str(location)+'\t1')
        else:
            pass


def process_trajectories(line):

    time_traj,lon_traj,lat_traj=get_traj_from_line(line)

    cutted_time,cutted_lon,cutted_lat=get_cutted_trajectory(time_traj,lon_traj,lat_traj)
           
    for cutted_time_traj,cutted_lon_traj,cutted_lat_traj in zip(cutted_time,cutted_lon,cutted_lat):

        refined_time_traj,refined_lon_traj,refined_lat_traj=get_refined_trajectory(cutted_time_traj,cutted_lon_traj,cutted_lat_traj)

        map_pop(refined_time_traj,refined_lon_traj,refined_lat_traj)

    
for line in read_from_text(sys.stdin):
    process_trajectories(line)
