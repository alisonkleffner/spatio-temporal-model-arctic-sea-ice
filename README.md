# A Spatio-Temporal Model for Arctic Sea Ice

### Alison Kleffner, Susan VanderPlas, and Yawen Guan

### University of Nebraska at Lincoln

<br>

Manuscript under preparation can be found [here](https://alisonkleffner.github.io/spatio-temporal-model-arctic-sea-ice/spatio-temporal-model-arctic-sea-ice.pdf)


**Abstract** 

Arctic Sea Ice is a barrier between the warm air of the ocean and the atmosphere, thus playing an important role in the climate. When narrow linear cracks (leads) form in the sea ice, the heat from the ocean is then released into the atmosphere. To estimate where cracks may form, motion data from the RADARSAT Geophysical Processing System (RGPS) are analyzed. The RGPS provides a set of trajectories (or cells) to trace the displacements of sea ice, however, chunks of data are missing due to the data collection method. We propose a spatial clustering and interpolation method that allows us to infer missing observations and estimate, where a crack may form. To do this feature inputs were created for KNN clustering by creating a bounding box around each trajectory, resulting in trajectories being assigned a cluster. A crack is considered to have formed on the boundary between different clusters. Within the clusters, spatiotemporal interpolation method is used to infer missing locations. Our clustering approach is then compared to other methods to determine ice crack formation, and cross-validation is used to assess our interpolation method.

