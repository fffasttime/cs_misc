create database stu;

create table info (name CHAR(6), age int);
insert into info values('fff',0);
insert into info values('dec',10);
insert into info values('oxa',16);
insert into info values('hex',16);
insert into info values('uuz',1000);
;;;


create table TeSt (tt iNt);
select * from __dbmeta
select * from __tbmeta_info
sElect * from __tbmeta_test
select * froM info where age<30
select * from info where name='fff' and 'fff'=name

select * from info where age>10
select * from info where age<=10
select * from info where age>10 and age<30
select * frOm INFO where age>10 or age<30
select * from info where name!='fff'
select * from info where name<>'uuz' and name!='fff'

select * from info WHERE age>10 or age<30 and name<>'hex'
select * from info WHERE (age>10 or age<30) and name<>'hex'


create table score (name CHAR(6), score int)
insert into score values('fff',0)
insert into score values('fff',1)
insert into score values('fff',2)
insert into score values('fff',3)
insert into score values('fff',4)
insert into score values('oxa',99)
insert into score values('oxa',100)

select * from info,score where score<3 and age<5

exit
