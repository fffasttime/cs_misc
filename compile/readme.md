# 简单的sql DBMS

## 简述

一个简化的关系数据库实现，同时也是编译原理课上机作业。

sql解析使用了lex(flex)+yacc(bison)。数据库部分使用c++编写。数据库运行时数据被读取到内存中，退出时可以保存到磁盘文件上。

支持select, insert, delete, update, drop, create等几个数据库常见操作。

数据的域可定义为 32位整数INT，定长字符串CHAR(N)，浮点数FLOAT 类型。在执行语句过程中，发现类型不匹配时会报错。

在功能上，还实现了多表连接查询，NULL类型，UNIQUE和NOT NULL完整性约束等。部分支持WHERE子句中的IN (select ...)嵌套查询功能。test/featrue.sql中给出了能正确运行示例语句。不过这些功能可能还有潜在的缺陷，而且执行效率有待优化。

## 说明

### 运行

环境：Linux (Ubuntu 16.04)，安装gcc, flex, bison

编译：make

若编译simplesql.cc出错，可尝试去掉-DGLIBCXX_DEBUG选项或更新gcc。

运行：

./simplesql  进入交互运行模式。

./simplesql  test/create.sql  会依次读取执行create.sql里的语句，若出错报错并继续执行下一条。

### 文件说明

simplesql.l 是lex文件，定义词法规则。

simplesql.y 是yacc文件，定义语法规则和语法动作。执行动作时一般会调用simplesql.h中的函数。

simplesql.h/simplesql.cc 包含对simplesql.y中用到的数据结构的定义。同时实现数据库的主要操作(select, insert, ...)。

storage.h/storage.cc 里包含了数据在内存中的格式的定义，以及如何从磁盘中读取和存储。

common.h/common.cc 包含少量公用的定义。

test/create.sql 包含了对于创建、查询、删除这几个基本功能的测试。

test/feature.sql 包含了几个测试中的功能，包括空值NULL，嵌套查询，UNIQUE约束等，均能正确运行。

data/* 数据文件保存在此处。

## 特性和功能

SQL语法对关键字和名称的大小写不敏感，所以引号外的字母大小写可随意。

针对输入命令中常见的错误，例如语法错误，未找到字段名，类型不匹配等，会给出错误提示，并忽略这条语句，继续执行下一条。

运行时过程中每行执行一条语句，可以不接分号。如果要在同一行执行多条语句，可以使用分号隔开。语句中/**/中和--后面的内容可被视为注释忽略。

### 创建、查看、保存数据库

数据库文件保存在./data目录下，文件名为”数据库名.db“。

一开始，命令提示符为sql>，此时需要创建一个新的数据库或加载已有数据。

#### 创建

使用CREATE DATABASE name; 或CREATE SCHEMA name;创建一个新的数据库。

执行这条语句后会创建一个名为__dbmeta的空表，这个表保存的是数据库内所有表的信息。

#### 加载

USE stu; 语句会读取./data/stu.db中的数据。此时命令提示符会变为stu>。

#### 查看

列出所有表的信息：

SHOW TABLES

列出某个表所有列的信息：

SHOW COLUMNS FROM table_name

由于在实现中，表的元信息也使用表保存，所以执行上述两条SHOW语句和使用如下select查询是等价的：

SELECT * FROM __dbmeta;  

SELECT * FROM __tbmeta\_table\_name;

#### 保存

注意由于数据操作在内存中进行，所有的更改操作不会立即写入磁盘。

使用save命令可保存当前内存中的数据。使用exit命令可以正常保存退出。

如果是Ctrl-C中断或者异常退出，则数据不会被保存。

### 查询

SELECT语句用于查询表的内容，例如：

查找表1和2连接后所有满足要求的列。

SELECT * FROM table1, table2 WHERE conditions

查找后选择两列。如果WHERE条件句或选择列中遇到同名的列，可以指定列所在的表名。查询结果的列名可以使用AS重命名。

SELECT table1.col1 AS 'new col name', table2.col2  FROM table1, table2 WHERE conditions

查询结果会在命令行中以表状显示，形如：

```
+----+----+-----+
|name|year|score|
+----+----+-----+
| fff|2020|   80|
| fff|2019|   75|
+----+----+-----+
found 2 items

```

SELECT实现上可支持多表连接查询，以及where子句中与外部无关的嵌套查询，但是效率不高。

where 子句中的二元运算符可为 AND, OR, =, >, <, >=, <=, !=, IN 。

其中，AND/OR左右只能连接比较的结果，比如 stu.name=score.name AND score.score>90 。

比较运算符 =, >, <, >=, <=, != 两侧可以都是列名称，也可以一侧是要比较的值。

IN用于嵌套查询，左侧是一个列名称，右侧是一个SELECT子查询。

### 数据修改

INSERT、DELETE和UPDATE用于数据修改。

INSERT INTO table VALUES (val1, val2, ...) 可插入一条数据，括号中的数据必须和表定义时的各个列相匹配。

INSERT INTO table (col1, col2, ...) VALUES (val1, val2, ...) 可指定要插入的列，未选择的列会被填充为NULL。

DELETE * FROM table; 或 DELETE FROM table; 会删除表中的所有数据，但保留当前表 。

DELETE FROM table WHERE condition; 删除表中满足condition条件的行。

UPDATE SET col_name=val WHERE condition; 将表中满足condition条件行中对应列的值修改为val。

### 表创建、删除

CREATE TABLE table_name(col1 CHAR(10) UNIQUE, col2 INT NOT NULL, col3 FLOAT);

可创建一个名为table_name的表，包含三个列col1, col2, col3。

第一列名称为col1，为定长10个字符的类型，而且要求唯一。第二列col2为NOT NULL的整数类型。

DROP TABLE table_name; 在当前模式中删除整个表。



## 部分要点

### 数据存储格式

数据存储结构在storage.cc/storage.h中实现。

一个database(schema)被保存为一个单独的.db文件，每个schema可以包含多个表table。

在创建或读取数据库时，一个名为__dbmeta的表首先被创建/读取，包含其他所有表的信息。

然后创建/读取每一个普通表时，会同时创建/读取有一个__tbmeta\_表名称 的表，记录它各个列的类型定义。

对于表的内容，按照二维表先行后列密集型方式存储。每一个域都有确定的长度，例如int占4个字节，定长字符串char(n)占固定的n+1个字节长度。即使空缺也会占用同样大小的空间。

这种结构缺点是相比现代数据库查询、修改的场景的效率不高，不过实现起来相对简单。

### 查询

查询过程主要在simplesql.cc/simplesql.h的Searcher类中。

一般来说，查询按照连接，选择，投影的步骤进行。

对于多个表的查询，Searcher在判断一个二元比较时会枚举所有使条件为真的行。准确的说，查询判断过程会枚举出满足条件的元组T，这些元组T表示几个表笛卡尔积的中满足条件的子集。由于未对查询语法树做代数优化，查询的中间条目也可能达到所有表的笛卡尔积，因此效率不高。

在查询到所有满足条件的条目后，再根据选择的列做投影，结果生成一个新的表。



## 声明

虽然test中的sql代码可正常执行，但是未做更多的测试，可能包含很多潜在的缺陷。因此上述代码仅用于课程学习的用途。

本代码按照MIT License开源，可以任意复制和修改。但作为课程作业不建议直接复制粘贴。

其中lex和yacc文件，语法树数据定义可作为参考。数据存储部分和查询部分由于设计原因可维护性较差，建议自己重新设计编写。