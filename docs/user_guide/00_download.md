title: Download

# How to download FEST-3D
The FEST-3D solver is available at [GitHub](https://github.com/FEST3D/FEST-3D) to download.  The tutorials are provided in separate [GitHub](https://github.com/FEST3D/run) repository and the same is used as submodule in main source code with folder name ``[RootFolder]/run``. The grid files for the tutorials are provided in the respective test case folders.

@note
In case you do not have Git installed on you local machine, you can do so by using following command on ubuntu: <br>
sudo apt-get install git<br>
A similar command can be used for other Linux packages.


### Use following to clone the git repository without tutorials
```
git clone https://github.com/FEST3D/FEST-3D.git
```

@note
The tutorials are submodule to main FEST-3D source code on the GitHub.

### Use following command to get git repository with tutorials
```
git clone --recursive https://github.com/FEST3D/FEST-3D.git
```
#### If you already have source code, then use following to downlaod tutorials
```
git submodule update --init
```

If you are downloading the FEST-3D package as Zip file from the [GitHub](https://github.com/FEST3D/FEST-3D) or [Home page](https://fest3d.github.io) on Linux folder use Unzip:
```
sudo apt-get install unzip
```
and use following command to unzip the folder
```
unzip FEST-3D-master.zip -d FEST3D
```
The run folder is not part of the main zip folder.
