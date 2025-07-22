# Converting IHP Open PDK for use with VACASK

Clone IHP SG13G2 PDK
```
git clone https://github.com/IHP-GmbH/IHP-Open-PDK
```
Now you have a directory named IHP-Open-PDK. 

You will need the path to VACASK's Python scripts. If you don't know where these scripts are, type
```
vacask -dp
```
and look for "Python path addition". Suppose the python path addition is "/usr/local/lib/vacask/python" and the PDK is in "/home/myname/IHP-Open-PDK". Type
```
IHPPDK=/home/myname/IHP-Open-PDK PYTHONPATH=/usr/local/lib/vacask/python python3 -m sg13g2tovc
```
The converter will process the Ngspice models and create a directory named "ihp-sg13g2/libs.tech/vacask" with the converted models. 





