# Incoherent Target Scattering Decomposition of Polarimetric SAR Data Based on Vector Model Roll-Invariant Parameters (Aghababaei Decomposition)


In this repository you can find the matlab code associated to the IEEE TGRS paper "Incoherent Target Scattering Decomposition of Polarimetric SAR Data Based on Vector Model Roll-Invariant Parameters", DOI: 10.1109/TGRS.2016.2540807 


This method extracts roll-invariant features from polarimetric SAR images. 
Its structure resembles Touzi decomposition, but with the below main distinction:
- Touzi decomposition extracts features that are roll-invariant to the orientation angle of the polarization ellipse of maximum co-polarization backscattering, which is descriptive of the polarization.
- Aghababaei decomposition extracts features that are roll-invariant to the orientation of the symmetry axis of the maximum symmetric component of the scatterer, which is descriptive of the target.

The code has meant for reaerch purpose If you use it, please cite as the following:


H. Aghababaee and M. R. Sahebi, "Incoherent Target Scattering Decomposition of Polarimetric SAR Data Based on Vector Model Roll-Invariant Parameters," in IEEE Transactions on Geoscience and Remote Sensing, vol. 54, no. 8, pp. 4392-4401, Aug. 2016, doi: 10.1109/TGRS.2016.2540807.


# Usage

You can find:
1) the matlab code of the method

Aghababaei_decomposion.m

2) A free available Flevoland testing image

    Feloveland.mat

3) A demo for running the code on the available data:

   Demo.m
