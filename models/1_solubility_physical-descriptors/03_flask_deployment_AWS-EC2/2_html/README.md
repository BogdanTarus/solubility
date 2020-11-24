
# instructions for linux/MacOS

## 1. Create an Ubuntu instance on EC2-AWS

## 2. Open two terminals on the local machine

## 3. Define in each terminal the privateKey, addressIP, and addressDNS:

privateKey="path/to/the/privatekey.pem"

> Examples of addressIP, and addressDNS. Change accordingly.

addressIP="3.88.33.133"
addressDNS="ec2-3-88-33-133.compute-1.amazonaws.com"


## 4. On the first terminal, ssh to EC2 instance :

ssh -i ${privateKey} ubuntu@${addressDNS}

## 5. On the second terminal, copy the files from the local machine to the EC2 instance

cd to/working/directory

scp -i ${privateKey} -r * ubuntu@${addressDNS}:/home/ubuntu

## 6. On the Ubuntu EC2 instance terminal, update and install requirements:

sudo apt-get update && sudo apt-get --yes install python3-pip
sudo pip3 install -r requirements.txt

sudo apt-get install python3-rdkit librdkit1 rdkit-data

## 7. launch it with

python3 app.py


## 8. Meanwhile, from another machine web browser, connect to the application by typing the adress (change the addressDNS accordingly)

ec2-3-88-33-133.compute-1.amazonaws.com:8080

Enter the SMILES structure of the compound (e.g. FC(F)(Cl)C(F)(Cl)Cl )

> The flask application will detect whether the SMILES syntax is correct.
