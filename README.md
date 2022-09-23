# A Spatio-Temporal Model for Arctic Sea Ice

### Alison Kleffner, Susan VanderPlas, and Yawen Guan

### University of Nebraska at Lincoln

<br>

Manuscript under preparation can be found [here](https://alisonkleffner.github.io/spatio-temporal-model-arctic-sea-ice/spatio-temporal-model-arctic-sea-ice.pdf)


**Abstract** 

Arctic Sea Ice is a barrier between the warm air of the ocean and the atmosphere, thus playing an important role in the climate. When narrow linear cracks (leads) form in the sea ice, the heat from the ocean is released into the atmosphere. To estimate where cracks may form, motion data from the RADARSAT Geophysical Processing System (RGPS) is analyzed. RGPS provides a set of trajectories of points on the ice sheet to trace the displacement of sea ice, however, chunks of data are missing due to the data collection method by satellite. We propose a spatio-temporal clustering and interpolation method that estimates where a crack may form and allows us to infer missing observations. Features based on the sea ice displacement are created for k-means clustering by creating a bounding box around each trajectory, resulting in trajectories being assigned a cluster. A crack is considered to have formed on the boundary between different clusters. Within the clusters, a spatio-temporal interpolation model using a Gaussian Process is used to infer missing locations. Our clustering approach is compared to ice lead detection methods, and cross-validation is used to assess our interpolation method.


