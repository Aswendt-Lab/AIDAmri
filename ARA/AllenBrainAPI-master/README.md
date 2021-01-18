# README #

Tools for accessing the Allen Atlas API from MATLAB. 
e.g. downloading images, projection data, searching for particular experiments, etc. 
Functionality isn't complete right now and is a bit rough in places, but everything works and the functions should provide a good template for showing you how to add functionality. 
For more details see the [Allen API documentation](http://help.brain-map.org/display/mouseconnectivity/API). 
You will need to install [JSONlab](http://ch.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files-in-matlab-octave) for some functions to work. 



### Examples ###

Get structure list that underlies the Allen Reference Atlas:
```
>> S=getAllenStructureList;
```


Keep only the areas that are children of `Cerebellum`:
```
>> S=getAllenStructureList('childrenOf','Cerebellum');
```

### Also see ###
There are other useful functions too. See:

* DownloadImageSeries.m 
* findAllenExperiments.m
* getInjectionIDfromExperiment.m
* getProjectionDataFromExperiment.m