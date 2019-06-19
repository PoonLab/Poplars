# Installing *Poplars*

* [Linux]()
* [Windows]()
* [Mac OS]()

##For Linux 
To install *Poplars* on a Linux system, we recommend cloning the source code from this repository. If you do not have superuser privileges, you may install *Poplars* from the [Python Package Index](https://pypi.org/) (`pip`).

###From Source 

1. Clone the GitHub repository 
    ```console
    kwade4@Jesry:~/git$ git clone https://github.com/PoonLab/Poplars.git
    Cloning into 'Poplars'...
    remote: Enumerating objects: 380, done.
    remote: Counting objects: 100% (380/380), done.
    remote: Compressing objects: 100% (226/226), done.
    remote: Total 380 (delta 191), reused 309 (delta 143), pack-reused 0
    Receiving objects: 100% (380/380), 30.65 MiB | 31.32 MiB/s, done.
    Resolving deltas: 100% (191/191), done.
    ```
    
    This will create a directory (folder) called `Poplars`
    
2. Navigate to the `Poplars` directory. 

    Then run the `setup.py` installation script with superuser `sudo` privileges, which will install the package and the system-specific MAFFT binaries (executables) to `/usr/local/lib`. This allows you to use *Poplars* anywhere you can use Python! 
    
    ```console
    kwade4@Jesry:~/git$ cd Poplars
    kwade4@Jesry:~/git/Poplars$ sudo python3 setup.py install
    [sudo] password for kwade4: 
    running install
    running build
    running build_py
    [...]
    Changing permissions of /usr/local/lib/python3.6/dist-packages/poplars/ref_genomes/K03455.fasta to 755
    Changing permissions of /usr/local/lib/python3.6/dist-packages/poplars/ref_genomes/M33262.fasta to 755
    Changing permissions of /usr/local/lib/python3.6/dist-packages/poplars/ref_genomes/HIV1_Mgroup.fasta to 755
    ```

###Using `pip`

1. To install *Poplars*, check which version of Python is on your system:
    
     ```console
        kwade4@Jesry:~$ python --version
        Python 2.7.15rc1
        kwade4@Jesry:~$ python3 --version
        Python 3.6.7
     ```
    
2. Ensure that the [Python Package Index](https://pypi.org/) is installed on your system. If you have Python 2, you will use `pip` or if you have Python 3, use `pip3`. 

    There are 2 ways to check if `pip` or `pip3` is installed:
    
    * Using the `which` command:
    
        If `pip`/`pip3` installed on your system, the location of the `pip`/`pip3` executable file(s) is displayed. 
        
        ```console
        kwade4@Jesry:~$ which pip3
        /usr/bin/pip3
        ```
        
        If `pip` or `pip3` is not installed on your system, nothing will be displayed this command will display nothing. 

        ```console
        kwade4@Jesry:~$ which pip3
        kwade4@Jesry:~$
        ```
    
    * Using the `pip -V` or `pip3 -V` command:
    
        If `pip` is installed on your system, the version of `pip` will be displayed.
    
        ```console
        kwade4@Jesry:~/git/Poplars$ pip3 -V 
        pip 9.0.1 from /usr/lib/python3/dist-packages (python 3.6) 
        ```
        
        If `pip` or `pip3` is not installed on your system, an error message will be printed. 
        
        ```console
        kwade4@Jesry:~$ pip3 -V
        pip3: command not found
        ```
        
    The exact output of these commands varies with the operating system. This example uses Ubuntu 18.0.4. 

3. If `pip` is not installed, you need to install `pip` using your system's package manager.     
    For example, installing `pip3` using Ubuntu 18.0.4: 
    ```console
    kwade4@Jesry:~$ sudo apt install python3-pip
    [sudo] password for kwade4: 
    Reading package lists... Done
    Building dependency tree       
    Reading state information... Done
    The following additional packages will be installed:
      dh-python libexpat1-dev libpython3-dev libpython3.6 libpython3.6-dev libpython3.6-minimal
      libpython3.6-stdlib python-pip-whl python3-dev python3-setuptools python3-wheel python3.6
      python3.6-dev python3.6-minimal
    Suggested packages:
      python-setuptools-doc python3.6-venv python3.6-doc binfmt-support
    The following NEW packages will be installed:
      dh-python libexpat1-dev libpython3-dev libpython3.6-dev python-pip-whl python3-dev python3-pip
      python3-setuptools python3-wheel python3.6-dev
    The following packages will be upgraded:
      libpython3.6 libpython3.6-minimal libpython3.6-stdlib python3.6 python3.6-minimal
    5 upgraded, 10 newly installed, 0 to remove and 189 not upgraded.
    Need to get 53.1 MB of archives.
    After this operation, 82.2 MB of additional disk space will be used.
    Do you want to continue? [Y/n] y
    ``` 
    
    Enter `y` if you are ready to install `pip3`.
    
    Note that this step requires a network connection. 
    
    ...
    

## For Windows 

To install *Poplars* on a Windows system, we recommend cloning the source code from this repository. 

###From Source 
1. Ensure [Git Bash](https://gitforwindows.org/) is installed. You can check this by searching for "Git Bash" in the Start Menu. 

2. Once Git Bash is installed, clone or download the repository. 

    


    
    
    
    





