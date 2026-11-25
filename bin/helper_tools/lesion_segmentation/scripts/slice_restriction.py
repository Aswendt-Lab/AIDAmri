"""Reference-free restriction to at most one 8-connected region per native z slice."""
import numpy as np
from scipy.ndimage import label


def one_region_per_slice(mask, seeds=()):
    """Choose a seed-hit region, otherwise largest; ties use deterministic label order.

    If seeds hit multiple regions, choose most seed hits, then largest area.
    Missed seeds never snap. Empty slices stay empty. This cannot separate an
    artifact joined to the lesion by foreground pixels.
    """
    mask = np.asarray(mask, dtype=bool)
    if mask.ndim != 3:
        raise ValueError('Expected a 3D mask in the native MRI grid')
    points = []
    for point in seeds:
        if len(point) != 3 or any(int(v) != v for v in point):
            raise ValueError('Seeds must be integer voxel coordinates')
        point = tuple(int(v) for v in point)
        if any(v < 0 or v >= mask.shape[k] for k, v in enumerate(point)):
            raise ValueError('Seed outside mask grid')
        if point not in points:
            points.append(point)
    result = np.zeros_like(mask)
    audit = []
    for z in range(mask.shape[2]):
        components, count = label(mask[:, :, z], structure=np.ones((3, 3)))
        areas = np.bincount(components.ravel(), minlength=count + 1)
        hits = np.zeros(count + 1, dtype=int)
        for i, j, k in points:
            if k == z and components[i, j]:
                hits[components[i, j]] += 1
        candidates = list(range(1, count + 1))
        chosen = max(candidates, key=lambda c: (hits[c], areas[c], -c)) if candidates else 0
        if chosen:
            result[:, :, z] = components == chosen
        audit.append(dict(slice=z,regions_before=count,regions_after=int(bool(chosen)),
                          selected_label=chosen,selected_by='positive_seed' if chosen and hits[chosen] else 'largest_area',
                          seeded_region_conflict=int(np.count_nonzero(hits)) > 1,
                          removed_voxels=int(mask[:, :, z].sum()-result[:, :, z].sum())))
    return result, audit


def continuous_seeded_region(mask, seeds=()):
    """One continuous slice-component path, anchored to the first positive hit.

    Neighbouring regions must touch under 26-connectivity: no blank-slice gaps
    or distant jumps. Dynamic programming favours seed hits, then path length,
    then total area. Other positive seeds constrain branch selection when
    reachable; missed or conflicting seeds are reported, never snapped.
    Without positives, the largest continuous path is used (unanchored).
    """
    mask = np.asarray(mask, dtype=bool)
    # Reuse coordinate validation; this does not select the final result.
    one_region_per_slice(mask, seeds)
    seeds = list(dict.fromkeys(tuple(int(v) for v in p) for p in seeds))
    planes=[];areas=[];votes=[]
    for z in range(mask.shape[2]):
        plane,n=label(mask[:,:,z],structure=np.ones((3,3)))
        planes.append(plane);areas.append(np.bincount(plane.ravel(),minlength=n+1))
        vote=np.zeros(n+1,dtype=int)
        for i,j,k in seeds:
            if k==z and plane[i,j]:vote[plane[i,j]]+=1
        votes.append(vote)
    edges=[]
    for z in range(len(planes)-1):
        a,b=planes[z:z+2];pairs=set();sx,sy=a.shape
        for dx in [-1,0,1]:
            for dy in [-1,0,1]:
                left=a[max(0,dx):min(sx,sx+dx),max(0,dy):min(sy,sy+dy)]
                right=b[max(0,-dx):min(sx,sx-dx),max(0,-dy):min(sy,sy-dy)]
                valid=(left>0)&(right>0)
                if valid.any():pairs.update(map(tuple,np.unique(np.stack([left[valid],right[valid]],axis=1),axis=0)))
        edges.append(pairs)
    hits=[int(planes[k][i,j]) for i,j,k in seeds]
    result=np.zeros_like(mask);chosen={}
    if seeds and not any(hits):
        return result,dict(seed_hits=[False]*len(seeds),warnings=['No positive seed intersects the network prediction. No lesion was selected; place another seed or review the network miss.'],selected_slices=[],anchor=None)
    def extend(z,c,direction):
        # Best anchored path ending in each reachable component at this layer.
        states={c:((int(votes[z][c]),1,int(areas[z][c])),[(z,c)])};best=next(iter(states.values()))
        k=z
        while 0<=k+direction<len(planes):
            nxt=k+direction;new={}
            for a,b in sorted(edges[min(k,nxt)]):
                parent,child=(a,b) if direction==1 else (b,a)
                if parent not in states:continue
                score,path=states[parent];candidate=((score[0]+int(votes[nxt][child]),score[1]+1,score[2]+int(areas[nxt][child])),path+[(nxt,int(child))])
                if child not in new or candidate[0]>new[child][0]:new[child]=candidate
            if not new:break
            states=new;k=nxt
            candidate=max(states.values(),key=lambda item:item[0])
            if candidate[0]>best[0]:best=candidate
        return best
    if seeds:
        index=next(i for i,c in enumerate(hits) if c);anchor=(seeds[index][2],hits[index])
        z,c=anchor
        for direction in [-1,1]:chosen.update(extend(z,c,direction)[1])
    else:
        # Deterministic largest continuous track across the volume.
        states={};best=((0,0),[])
        for z in range(len(planes)):
            new={c:((int(areas[z][c]),1),[(z,c)]) for c in range(1,len(areas[z]))}
            if z:
                for parent,child in sorted(edges[z-1]):
                    if parent not in states:continue
                    score,path=states[parent];candidate=((score[0]+int(areas[z][child]),score[1]+1),path+[(z,int(child))])
                    if candidate[0]>new[child][0]:new[child]=candidate
            states=new
            if states:
                candidate=max(states.values(),key=lambda item:item[0])
                if candidate[0]>best[0]:best=candidate
        chosen.update(best[1]);anchor=None
    for z,c in chosen.items():result[:,:,z]=planes[z]==c
    retained=[bool(result[p]) for p in seeds];warnings=[]
    if not seeds:warnings.append('Avoidance-only selection is unanchored: the largest continuous track may be the wrong region. A positive seed is recommended.')
    if seeds and not all(retained):warnings.append('Some positive seeds were missed or cannot belong to one continuous region per slice. Review the seed-hit indicators; add avoidance ROIs or revise conflicting seeds.')
    return result,dict(seed_hits=retained,warnings=warnings,selected_slices=sorted(chosen),anchor=anchor,
                       rule='one 8-connected region per native z slice; adjacent regions must be 26-connected; no gaps or snapping')
