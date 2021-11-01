#!/usr/bin/env python3

## Create CSV file to be read by pandas?
## Easy to use data science/machine learning machinery on it

import mysql.connector

mydb = mysql.connector.connect(
  host="sql3.freesqldatabase.com",
  user="sql3448074",
  password="JWiRYDR581",
  port="3306"
)

mycursor = mydb.cursor()

