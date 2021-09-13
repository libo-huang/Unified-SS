# Unified-SS

## A Unified Optimization Model of Feature Extraction and Clustering for Spike Sorting

The MATLAB implementation of automatic feature extraction and clustering algorithm for spike sorting, as proposed in **Huang, Libo, Lu Gan, and Bingo Wing-Kuen Ling. "A Unified Optimization Model of Feature Extraction and Clustering for Spike Sorting." IEEE Transactions on Neural Systems and Rehabilitation Engineering 29 (2021): 750-759.** (available at [https://ieeexplore.ieee.org/abstract/document/9408665/](https://ieeexplore.ieee.org/abstract/document/9408665/))

If you find this algorithm useful in your work, please give credit by citing the above paper. If you have any questions regarding the algorithm or implementation, please do not hesitate to write to the author with the following email address: www(dot)huanglibo(AT)gmail(dot)com

You need MATLAB software to use this program.

## Usage
You can easily run the `demo` of the Unified-SS with (the 'current folder' path is Unified-SS),
```Matlab
>>> demo
```
in which two parts of the test are included, `fix k=3` and `estimate class number with 4 methods`. 

Specifically, the `fix k=3` part means the demo codes are implemented with fixed class/neurons number with 3 followed by the other spike sorting processes including spike filtering, spike detection, and the proposed unified feature extraction and clustering.

On the other hand, the `estimate class number with 4 methods` part has tested an automatic spike sorting framework with estimating the class/neurons number and unified feature extraction and clustering, etc.

For more details, please refer to the paper **Huang, Libo, Lu Gan, and Bingo Wing-Kuen Ling. "A Unified Optimization Model of Feature Extraction and Clustering for Spike Sorting." IEEE Transactions on Neural Systems and Rehabilitation Engineering 29 (2021): 750-759.** (available at [https://ieeexplore.ieee.org/abstract/document/9408665/](https://ieeexplore.ieee.org/abstract/document/9408665/)).

## Citation
```
@article{huang2021unified,
  title={A Unified Optimization Model of Feature Extraction and Clustering for Spike Sorting},
  author={Huang, Libo and Gan, Lu and Ling, Bingo Wing-Kuen},
  journal={IEEE Transactions on Neural Systems and Rehabilitation Engineering},
  volume={29},
  pages={750--759},
  year={2021},
  publisher={IEEE}
}
```

## Licence
Copyright(C) 2021, Libo Huang, www.huanglibo@gmail.com

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
