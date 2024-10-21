import numpy as np


from sklearn.cluster import DBSCAN



def dbscan_f(df, epsilon, minpts):


    X = np.vstack((df['x'],df['y'])).T
    
    # perform dbscan
    db = DBSCAN(eps=epsilon, min_samples=minpts).fit(X)
    group = np.int32(db.labels_)

    # save cluster
    df['cluster_ID'] = group
    df_cluster = df.loc[df['cluster_ID'] != -1]
    
    df_cluster = df_cluster.sort_values(by = ['cluster_ID'])

    '''
    Generating rec array file for picasso render with all localizations assigned to a cluster
    Colorcoding = protein_ID
    
    '''
    
    LOCS_DTYPE = [
        ('frame', 'u4'),
        ('x', 'f4'),
        ('y', 'f4'),
        ('photons', 'f4'),
        ('sx', 'f4'),
        ('sy', 'f4'),
        ('bg', 'f4'),
        ('lpx', 'f4'),
        ('lpy', 'f4'),
        ('group', 'u4'),
        # ('cluster_ID', 'u4')
    ]
    db_locs_protein_ID = np.rec.array(
        (df_cluster.frame, df_cluster.x, df_cluster.y, df_cluster.photons, df_cluster.sx, df_cluster.sy, df_cluster.bg,
         df_cluster.lpx, df_cluster.lpy, df_cluster.cluster_ID), dtype=LOCS_DTYPE,
    )
    
    
    return db_locs_protein_ID


def dbscan3d_f(df, epsilon, minpts):

    pxsize = 130 # nm
    X = np.vstack((df['x'],df['y'],df['z']/pxsize)).T
    
    # perform dbscan
    db = DBSCAN(eps=epsilon, min_samples=minpts).fit(X)
    group = np.int32(db.labels_)

    # save cluster
    df['cluster_ID'] = group
    df_cluster = df.loc[df['cluster_ID'] != -1]
    
    df_cluster = df_cluster.sort_values(by = ['cluster_ID'])

    '''
    Generating rec array file for picasso render with all localizations assigned to a cluster
    Colorcoding = protein_ID
    
    '''
    
    LOCS_DTYPE = [
        ('frame', 'u4'),
        ('x', 'f4'),
        ('y', 'f4'),
        ('z', 'f4'),
        ('photons', 'f4'),
        ('sx', 'f4'),
        ('sy', 'f4'),
        ('bg', 'f4'),
        ('lpx', 'f4'),
        ('lpy', 'f4'),
        ('group', 'u4'),
        ('cluster_ID', 'u4')
    ]
    db_locs_protein_ID = np.rec.array(
        (df_cluster.frame, df_cluster.x, df_cluster.y, df_cluster.z, df_cluster.photons, df_cluster.sx, df_cluster.sy, df_cluster.bg,
         df_cluster.lpx, df_cluster.lpy, df_cluster.protein, df_cluster.cluster_ID), dtype=LOCS_DTYPE,
    )
    
    
    return db_locs_protein_ID



