create schema feature

create table tb(name CHAR(5) UNIQUE, val FLOAT NOT NULL, val2 INT)

show databases;
show tables;
show columns from aaa; -- table not found
show columns from tb;

insert into tb values(NULL, NULL, 1) -- error: name can't be null
insert into tb values('s1', NULL, 1) -- error: val can't be null
insert into tb values('ssssss', 3.5, 1) -- error: string too long
insert into tb values('sssss', 3.5, 1) 
insert into tb values('s1', 3.5, NULL) 
insert into tb values('/**/', 3.5, 2) 
insert into tb values('--f', 4, 3) 

select * from tb 
select * from tb where val=3.5 -- 3 items
select * from tb where val=4 -- 1 item
select * from tb where val<4 -- 3 items

insert into tb (val2, val) values(4, 3) -- error: name can't be null
insert into tb (name) values('s1') -- error: name should be unique
insert into tb (name, val) values('s1', 10) -- error: name should be unique

save
