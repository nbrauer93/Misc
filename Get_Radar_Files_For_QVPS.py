
### **** Source code is courtesy of Randy Chase *****

import s3fs
import numpy as np
import tqdm
import pandas as pd

def makedtimes(in_list,bucket,date,radar):
  """This is a helper function to help select times"""
  s = pd.Series(in_list)
  split_s = s.str.split(bucket + date.strftime('/%Y/%m/%d/') + radar + '/' + radar+date.strftime('%Y%m%d_'),expand=True)
  split_s = split_s[1]
  counts = split_s.str.len()
  if date < np.datetime64('2007-12-31'):
    version = '.gz'
    min_count = 9.0
  elif date < np.datetime64('2015-12-31'):
    version = '_V06'
    min_count = 13
  else:
    min_count=10.0
    if radar[0]=='K':
      version = '_V06'
    elif radar[0] =='T':
      version = '_V08'
  split_s = split_s.str.split(version,expand=True)
  split_s = split_s[0]
  dtime = pd.to_datetime(date.strftime('%Y-%m-%d ')+ split_s).values
  df = pd.DataFrame(in_list,index=dtime)

  df['counts'] = pd.Series(counts.values,index=dtime)
  df = df.where(df.counts == min_count).dropna()
  df = df.drop('counts',axis=1)
  return df

#this is where the radar data are held
bucket = 'noaa-nexrad-level2'
#pick your favorite date
date = pd.to_datetime('2022-09-28')
#pick your favorite radar
radar_name = 'KTBW'
# Use the anonymous credentials to access public data
fs = s3fs.S3FileSystem(anon=True)

#Grab a list of all files on this day
in_list = fs.ls('s3://'+bucket + date.strftime('/%Y/%m/%d/') + radar_name)

#use our defined function to add in datetimes to allow for easy slicing
df = makedtimes(in_list,bucket,date,radar_name)

#select start hour UTC
s_hour = '1200'
#select end hour UTC
e_hour = '2359'

s_time = date.strftime('%Y-%m-%d ') + s_hour
e_time = date.strftime('%Y-%m-%d ') + e_hour
df_wanted = df.loc[slice(s_time,e_time)]
list_o_files = df_wanted.values



# download
savedir = './20220928/'
for file in tqdm.tqdm(list_o_files):
    fs.get(file[0], savedir+file[0].split('/')[-1])

#grab the filenames
import glob
rad_files = glob.glob(savedir + radar_name+'*')
rad_files.sort()
