# internal mode  
Storage mode should detach from logic mode when design, for convinience to modify its structure.  
## disk storage  
Here the data table is simply united by lines(aka records). It is not so efficiency now, but esaier.  
A more efficient way for query is save it by columns. A modern way for bigdata is to save key-value in hashmap, see Google Bigtable.  
We need an ADT support inserting and deleting operation on records. Specifically, fast insertion, deleting and enumerating.  
Data strcuture like B-tree is appropate.  
To simplify the implemention here, the data table will be fully loaded to memory at work, and save to disk when exit.  
So a dynamtic array (c++ std::vector) is enough.  
## table storage structure  
Before every data table(inculde metadata), a tid number indicates which table follow.  
tid is unique, and database metadata always have tid=0.  
The length and type of each field is loaded in metadata before, so the full length is known(dbmetadata.items*tbmetadta.fieldlength)  
## data field  
All data fields have fixed length, so it's easy to determine the position.  
(Todo) Save VARCHAR or BLOB object as a pointer in main table.  

# metadata  
Metadata saves in same type tables, similar as ordinary data.  
Metadata table always get full loaded in memory.  
## database metadata  
Each database have one table, always occupy the head of .db file.  
First 64 Byte is metadata of full database  
|DB_identification("SIMPLEDB") 10 byte | version 10 byte | other infomathon(reserved)| number of tables|  
Then a table lists ordinary tables in this database    
{'name':char(32), 'count':int, 'count_field':int}  
## table matadata  
Then every table have a table save its field info  
{'type':int, 'extra':int, 'name', (TODO)'foreign_key'}  

# selection
Here selection is simply enumerate all possible combinions, it's easy but very slow.  
Calcutations between tables clause " select * from some_table where a + b < c " can't be represented under this way.

# test & debug  
sql interpreter should be tested by examples.  
