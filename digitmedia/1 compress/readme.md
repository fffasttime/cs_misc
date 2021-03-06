## usage  
ArithmeticCoding.exe encode \<filename\>  
or  
ArithmeticCoding.exe decode \<filename\>  



## 算术编码

### 算法描述

算术编码是现在常用的压缩编码方法，它比霍夫曼编码更优，更接近信息熵。

把整个信源序列表示成0-1内的一个区间，其长度等于该序列的概率，在区间内选择一个最小位数的小数作为编码。序列中每个新增元素会缩短区间，就需要更多位来表示。可以发现这种编码中符号的平均码长可以为小数，而哈夫曼编码必须为整数。

在不考虑精度问题时，这个过程是相当简单的。我们只需要用单个字符来更新区间下限L和区间长度R：

L’ = L + R*(S[i]/S) R’=R*(S[i+1]/S – S[i]/S)

这里S[i]/S即单个符号i所划分区间的起始位置。

### 规整化和再归一化

#### 规整化

为了在有限精度的机器上实现，需要对编码区间取近似值。例如，当使用8位编码，符号abc等概率时，划分为[0, 82/256), [82/256, 164/256), [164/256, 1)。

但划分几次后，这个过程就不能继续进行了。虽然可以使用高精度计算把划分的区间变得很细，但它的计算效率没有应用的意义。

注意到，当区间收窄时，二进制小数的最高位可以逐步确定。这样可以把高位作为编码结果输出，然后按位左移当前数值，就可以在固定字长的寄存器上计算了。这个过程需要再归一化。

#### 再归一化

然而有时编码的最高位一直难以确定。例如上述输入为bbbbbb时，区间一直包括0.5，这样仍会超出机器精度限制。

为此，我们不得不设置一个阈值，当区间范围R小于阈值时，强行将编码的高位输出。这样做的后果是，如果区间真的在0.5上方，左移之后会大于1，也称为进位。因此，我们需要记录下边界L中连续0.01111111…中1的出现长度，如果出现进位，则全部改用0.10000000…填充。

下面代码使用16位整数作为上限，将连续8位一起处理，以提高效率：

如果在更新阶段，发生L>>8 > 0xff时，说明发生溢出，此时标记z=1，buf++，然后删除溢出位。

然后进行重归一化：当剩余区间长度小于1<<8时，我们把尝试把高8位弹出。当高8位为0xff(11111111B)时，说明暂时不可确定，计数器+1。否则先输出进位前面buffer中的内容，然后根据是否进位输出连续的1或者0。

#### 解码

解码过程就相对简单了，我们用一个V表示当前读到的二进制流。如果V小于区间下界L，则说明此处编码时发生了上溢，要把1<<p加回去。其他部分代码和编码时基本相同。

### 测试

对一个116704byte (约120kb)大小左右的测试文件进行压缩。测试文件名为test.cpp，是一个c++程序代码，使用ASCII编码。上面的一列数字为编码的划分方案，会被单独存在一个字典文件中，解压时会被使用。

经计算，该文件的香农熵为78982.1 byte。根据理论，算术编码的压缩后大小为其香农熵的上取整。这里算术编码压缩后大小为78983byte，与之相符。



## LZW编码

字典压缩方法，部分考虑了字之间的相关性。

初始化字典集合为字符集。

设当前处理到的区间为P，新读入一个字符为C。

若PC未在字典中，则按字典输出P，将PC加入字典，再令P=C。

若PC已在字典中，则P=PC。

解压流程与之类似，读入编码后的数据，按流方式在构建字典的同时解压。这说明不需要保存字典。

代码实现中，字典项的代码逐步+1，当一字节不够时表示成两字节，再不够时表示成三字节这样。字节的最高位用来区分这是几字节的编码。

### 测试

116704byte (约120kb)大小左右的测试文件，压缩结果为68138byte。

对同样的文件进行LZW字典编码，压缩率高于算术编码。这个结果合理的，因为算术编码是一阶熵编码，只统计了频率而没有考虑到相关性。