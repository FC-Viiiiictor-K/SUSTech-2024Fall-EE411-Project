# EE411 信息论与编码 课程项目

## 项目总览

`50-SF.txt`是课程提供的已经经过预处理的DNA序列读取结果。

`50-SF-decoded.jpg`是经过解码之后得到的《神奈川冲浪里》图片。

`robust_soliton.py`来自课程提供的[示例代码](https://github.com/jeter1112/dna-fountain-simplified/blob/main/utils/robust_soliton.py)，其实现了稳定孤波分布（Robust Soliton Distribution，RSD）上的采样器。
- 稳定孤波分布用于确定LT码中各个码元的“度数”，即码元应为多少个源文件片段的异或和。
- 由于编码工作是在Python2中完成的，而`random.sample`函数在Python2和Python3中的实现并不相同，原始代码无法直接在Python3中用于解码。为了在Python3使用该代码，我们在其中增加了一个符合Python2中的`random.sample`函数行为的随机采样函数实现。

`main.py`是项目的主程序，为解码的主要代码。

## 运行项目

### 环境要求

- Python ~= 3.10.8

### 运行代码

首先运行该命令安装依赖：

```bash
pip install -r requirements.txt
```

然后运行主程序：

```bash
python main.py
```

运行主程序之前可以先删除`50-SF-decoded.jpg`，以展示生成的解码结果。

解码得到的图片为：![50-SF-decoded.jpg](50-SF-decoded.jpg)

运行结束后，主程序将会给出为解码图片而接收的**可恢复的**码元数量：

```bash
Full file is recovered after receiving 1821 droplets.
```

## 原理理解

### Reed-Solomon纠错码

#### RS码的基本思想与编码过程

RS码的本质是基于原始信息构造**有唯一解的超定方程组**。这种超定方程组是“冗余”的，也就是说，即使丢失了一些方程，我们也可以根据剩下的方程组求得唯一解。RS码最自然的用途是**纠删码**，也就是说，当码字上的一小部分项丢失的时候，我们可以用RS码将丢失的项全部恢复。在现实中，RS码被广泛应用于CD和DVD光盘、条形码和二维码、RAID6标准等场景。除了纠删码以外，经过精细的算法构造，RS码也可以胜任**纠错码**的角色，也就是在不了解错误位置的情况下，检测接收到的码字是否存在错误，以及对码字中出现错误的项进行修正。

RS码的编码通常是在有限域 $GF(2^8)$ 上进行的，因为计算机使用的数据单位通常都是byte，其有 $2^8$ 种可能的取值。在本项目中，RS码选择的有限域也是 $GF(2^8)$ 。没有选择“碱基所在的” $GF(2^2)$ 的一大原因是域的大小限制了多项式的次数，进而限制了最多能编码的信息长度及提供纠错能力的冗余编码的长度。

设我们希望编码的信息为 $M$ ，其长度为 $k$ ，也就是我们希望对 $\{M_0, M_1, M_2, ..., M_{k-1}\}$ 进行编码；同时，设冗余的纠错码长度为 $m$ ，那么接收端最多能够容忍丢失 $m$ 个项（RS码作为纠删码），或者有 $\lfloor\frac{m}{2}\rfloor$ 个项出现错误（RS码作为纠错码）。编码完成后，我们应当获得一个长度为 $n=k+m$ 的码字。设 $GF(2^8)$ 上的一个生成元为 $\alpha$ ，则有 $m$ 次生成多项式 $g(x)=\prod_{i=0}^{m-1}(x-\alpha^i)$ ；同时，定义 $k-1$ 次信息多项式 $M(x)=\sum_{i=0}^{k-1}M_ix^i$ ，令 $s(x)=M(x)x^m-(M(x)x^m\mod g(x))$ ，那么 $n-1$ 次多项式 $s(x)$ 的系数列 $\{s_0, s_1, ..., s_{n-1}\}$ 就是编码的结果。

- 不难发现， $M(x)x^m$ 的最低次项为 $x^m$ ，而 $M(x)x^m\mod g(x)$ ，也就是 $M(x)x^m$ 除以 $g(x)$ 的余数多项式，其最高次项必然比 $g(x)$ 的次数少 $1$ ，即  $m-1$ ；设这个余数多项式为 $r(x)=-(M(x)x^m\mod g(x))$ ，其系数为 $\{r_0, r_1, ..., r_{m-1}\}$ ，那么 $M(x)x^m$ 和 $r(x)$ 的各项系数就不会相互影响， $s(x)=M(x)x^m+r(x)$ 的系数列即为 $\{s_0, s_1, ..., s_{n-1}\}=\{r_0, r_1, ..., r_{m-1}, M_0, M_1, ..., M_{k-1}\}$ 。

为了更好的直觉理解，我们尝试用线性代数语言来描述编码过程。我们希望构造 $n\times k$ 的分块矩阵 $A'$ 使得：

$$
A'\vec{M}=
\begin{bmatrix}
A\\
I_k
\end{bmatrix}
\vec{M}
=\vec{s}=\begin{bmatrix}
\vec{r}\\
\vec{M}
\end{bmatrix}
=\begin{bmatrix}
r_0\\
r_1\\
\vdots\\
r_{m-1}\\
M_0\\
M_1\\
\vdots\\
M_{k-1}
\end{bmatrix}
$$

其中， $A$ 是一个 $m\times k$ 的子矩阵，其含义为“从 $M$ 生成余数多项式 $r$ ”这一过程；这样，我们就构造出了一个关于 $\vec{M}$ 的超定方程组。由于解该方程组时，接收端仅了解接收到的不完整的 $\vec{s}$ 以及算法中选择的 $\alpha$ ，**子矩阵 $A$ 应当仅依赖于 $\alpha$**。

显然 $s(x)$ 可被 $g(x)$ 整除，即存在 $k-1$ 次多项式 $q(x)$ 使得 $g(x)q(x)=s(x)$ ；设 $g(x)$ 的系数列为 $\{g_0, g_1, ..., g_m\}$ ， $q(x)$ 的系数列为 $\{q_0, q_1, ..., q_{k-1}\}$ ，则有：

$$
G\vec{q}=\begin{bmatrix}
g_0 & 0 & 0 & \dots & 0 & 0\\
g_1 & g_0 & 0 & \dots & 0 & 0\\
g_2 & g_1 & g_0 & \dots & 0 & 0\\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots\\
0 & 0 & 0 & \dots & g_m & g_{m-1}\\
0 & 0 & 0 & \dots & 0 & g_m\\
\end{bmatrix}
\begin{bmatrix}
q_0\\
q_1\\
\vdots\\
q_{k-1}
\end{bmatrix}
=\begin{bmatrix}
s_0\\
s_1\\
\vdots\\
s_{n-1}
\end{bmatrix}
=\begin{bmatrix}
r_0\\
r_1\\
\vdots\\
r_{m-1}\\
M_0\\
M_1\\
\vdots\\
M_{k-1}
\end{bmatrix}
$$

其中 $G$ 是一个 $n\times k$ 的矩阵，其含义为“乘以多项式 $g(x)$ ”；其构造方法为从左上角开始的对角线上全部为 $q_0$ ， $q_0$ 的下方全部为 $q_1$ ，以此类推，直到 $q_{k-1}$ 恰好填满结束于右下角的对角线。进一步地，我们将 $G$ 分解为两个块矩阵 $G_r$ 和 $G_M$ ：

$$
G=\begin{bmatrix}
G_r\\
G_M
\end{bmatrix}
$$

其中， $G_r$ 是一个 $m\times k$ 的矩阵， $G_M$ 则是一个 $k\times k$ 的**上三角矩阵**，那么有：

$$
G_M\vec{q}=\vec{M}\rightarrow \vec{q}={G_M}^{-1}\vec{M}\\
\vec{r}=G_r\vec{q}=G_{r}{G_M}^{-1}\vec{M}
$$

那么我们成功构造得到 $A=G_{r}{G_M}^{-1}$ 。

#### RS码的解码过程

在已知丢失项的位置或出错项的位置的情况下，纠错的算法是显然的；将RS码作为纠错码的重点是**如何找到出错的项的位置**。

设接收端收到的码字为 $R=\{R_0, R_1, ..., R_{n-1}\}$ ，对应的码字多项式为 $R(x)=\sum_{i=0}^{n-1}R_ix^i$ ；设接收端收到的码字和原信息的差为 $E=\{E_0, E_1, ..., E_{n-1}\}$ ，则有“错误多项式” $E(x)=R(x)-s(x)=\sum_{i=0}^{n-1}(R_i-M_i)x^i=\sum_{i=0}^{n-1}E_ix^i$ 。通过求得 $E(x)$ ，接收端可以计算 $s(x)=R(x)+E(x)$ ，由此实现解码。

由于 $g(x)q(x)=s(x)$ ，对于 $\forall i\in[0,m)$ ，将生成元的幂次 $\alpha^i$ 代入，则有 $s(\alpha^i)=g(\alpha^i)q(\alpha^i)=(\alpha^i-\alpha^0)(\alpha^i-\alpha^1)...(\alpha^i-\alpha^i)...(\alpha^i-\alpha^m)q(\alpha^i)=0$ ；那么将 $\alpha^i$ 代入 $R(x)$ ，则有 $R(\alpha^i)=E(\alpha^i)+s(\alpha^i)=E(\alpha^i)$ ，不妨记 $S_i=R(\alpha^i)=E(\alpha^i)$ ，则：

$$
\begin{bmatrix}
(\alpha^0)^0 & (\alpha^0)^1 & \dots & (\alpha^0)^{n-1}\\
(\alpha^1)^0 & (\alpha^1)^1 & \dots & (\alpha^1)^{n-1}\\
\vdots & \vdots & \ddots & \vdots \\
(\alpha^{m-1})^0 & (\alpha^{m-1})^1 & \dots & (\alpha^{m-1})^{n-1}\\
\end{bmatrix}
\begin{bmatrix}
E_0\\
E_1\\
\vdots\\
E_{n-1}
\end{bmatrix}
=\begin{bmatrix}
S_0\\
S_1\\
\vdots\\
S_{m-1}
\end{bmatrix}
$$

然而，该方程组**欠定**，无法求解；为了求解，接收端必须找到出错的项的位置。

设有 $\nu\in(0,\lfloor\frac{m}{2}\rfloor]$ 个位置出现错误，分别为 $I=\{I_0,I_1,...,I_{\nu-1}\}$ ，那么必然只有 $E_{I_0}$ 、 $E_{I_1}$ 、…、 $E_{I_{\nu-1}}$ 的值非零， $E$ 的其他系数全部为 $0$ 。不难发现，此时 $S_i=\sum_{j=0}^{\nu-1}E_{I_j}(\alpha^i)^{I_j}=E_{I_0}(\alpha^i)^{I_0}+E_{I_1}(\alpha^i)^{I_1}+...+E_{I_{\nu-1}}(\alpha^i)^{I_{\nu-1}}$ 。接下来，设**定位多项式** $\Lambda(x)=\prod_{i=0}^{\nu-1}(1-\alpha^{I_i}x)=1+\Lambda_1x+\Lambda_2x^2+...+\Lambda_\nu x^\nu$ ，则对于 $\forall i,j\in[0,\nu)$ 有：

$$
\begin{aligned}
\Lambda(\alpha^{-I_i})=&(1-\alpha^{I_0}\alpha^{-I_i})(1-\alpha^{I_1}\alpha^{-I_i})...(1-\alpha^{I_i}\alpha^{-I_i})...(1-\alpha^{-I_{\nu-1}}\alpha^{-I_i})\\
=&(1-\alpha^{I_0}\alpha^{-I_i})(1-\alpha^{I_1}\alpha^{-I_i})...(1-1)...(1-\alpha^{-I_{\nu-1}}\alpha^{-I_i})\\
=&0\\
E_{I_i}(\alpha^{I_i})^{j+\nu}\Lambda(\alpha^{-I_i})=&E_{I_i}(\alpha^{I_i})^{j+\nu}+\Lambda_1E_{I_i}(\alpha^{I_i})^{j+\nu-1}+\Lambda_2E_{I_i}(\alpha^{I_i})^{j+\nu-2}+...+\Lambda_{\nu}E_{I_i}(\alpha^{I_i})^{j}\\
=&0\\
\sum_{i=0}^{\nu-1}E_{I_i}(\alpha^{I_i})^{j+\nu}\Lambda(\alpha^{-I_i})=&\sum_{i=0}^{\nu-1}(E_{I_i}(\alpha^{I_i})^{j+\nu}+\Lambda_1E_{I_i}(\alpha^{I_i})^{j+\nu-1}+\Lambda_2E_{I_i}(\alpha^{I_i})^{j+\nu-2}+...+\Lambda_{\nu}E_{I_i}(\alpha^{I_i})^{j})\\
=&\sum_{i=0}^{\nu-1}E_{I_i}(\alpha^{I_i})^{j+\nu}+\Lambda_1\sum_{i=0}^{\nu-1}E_{I_i}(\alpha^{I_i})^{j+\nu-1}+...+\Lambda_{\nu}\sum_{i=0}^{\nu-1}E_{I_i}(\alpha^{I_i})^j\\
=&S_{j+\nu}+\Lambda_1S_{j+\nu-1}+...+\Lambda_{\nu}S_j\\
-S_{j+\nu}=&\Lambda_1S_{j+\nu-1}+...+\Lambda_{\nu}S_j
\end{aligned}
$$

那么，对于 $\forall j\in[0,\nu)$ ，有：

$$
\begin{bmatrix}
S_{\nu-1} & S_{\nu-2} & \dots & S_0\\
S_\nu & S_{\nu-1} & \dots & S_1\\
\vdots & \vdots & \ddots & \vdots\\
S_{2\nu-2} & S_{2\nu-3} & \dots & S_{\nu-1}
\end{bmatrix}
\begin{bmatrix}
\Lambda_1\\
\Lambda_2\\
\vdots\\
\Lambda_\nu
\end{bmatrix}
=\begin{bmatrix}
-S_\nu\\
-S_{\nu+1}\\
\vdots\\
-S_{2\nu-1}\\
\end{bmatrix}
$$

此时，**可以解得 $\vec{\Lambda}$**；之后，对于 $\forall i\in[0,n)$ ，将 $\alpha^i$ 代入生成多项式，若 $\Lambda(\alpha^i)=0$ ，则说明 $\exists j\in[0,\nu)$ 使得 $I_j=i$ 。这样，所有出错的项的位置就都被找到了。


### Luby Transform编码

#### 喷泉码

Luby Transform编码（Luby Transform Code，下称“LT码”），是世界上第一种基于喷泉码（Fountain Code）的编码方式。喷泉码是一种源编码技术，其设计目标是使得即使通信信道的噪声大、丢包率高，接收端也可以非交互地通过接收一段码元流来恢复出原始数据。在喷泉码中，发送端“像喷泉喷水一样”生成一个无限长的码元流，接收端通过“像接喷泉的水滴一样”接收足够多的码元，就可以恢复出原始数据。喷泉码的每个码元都是输入数据块的某些线性组合，接收者无需知道哪些码元是“关键”的，也无需专门请求特定码元，只要接收到足够数量的码元，就可以恢复原始数据。喷泉码的设计特点使得其特别适合“以数据包形式读取信息”、“丢包率较高”的场景，例如移动通信、卫星通信等，以及**DNA存储与读取**。

喷泉码与DNA读取的特点相符。在DNA读取的过程中，读到不完整或具有增删错误的DNA序列是不可避免且常见的（也就是“丢包率高”），而且无法指定读取特定的DNA序列。这非常符合喷泉码的设计目标——喷泉码的非交互性和丰富的冗余度使得其能够仅通过读取大量DNA序列来恢复原始数据。

#### LT码

在LT码中，发送端首先将原始数据分割成若干个数据块（Chunk），然后源源不断地基于这些数据块生成码元（Symbol，或者说“水滴/Droplet”）流，每个码元都包含一个或多个数据块的异或和（也就是“线性组合”，这个异或和或线性组合被称为“负载/Payload”）。除了包含数据块的异或和之外，每个码元还需要记录其包含了多少数据块的异或和（也就是这个码元的“度数”），以及是哪些数据块的异或和，以便接收端基于这些信息恢复所有数据块。

形式化地，设原始数据包含了 $n$ 个数据块 $\{c_1, c_2, ..., c_n\}$ ，发送端生成了 $m$ 个码元 $\{s_1, s_2, ..., s_m\}$ ；令 $m\times n$ 的01矩阵 $A$ 表示码元包含了哪些数据块的异或和，即 $A_{ij}=1$ 当且仅当 $s_i=...\oplus c_j\oplus...$ ；那么有：

$$
A\vec{c}=\vec{s}
$$

在LT码中，码元的度数是遵循稳定孤波分布的，也就是说，许多码元仅是相当少量数据块的异或和；那么，矩阵$A$通常是一个**稀疏矩阵**。

假设接收端最终接收到了 $m'$ 个码元 $\{s_{i_1}, s_{i_2}, ..., s_{i_{m'}}\}$ ，令：

$$
A' = \begin{bmatrix}
A_{i_1}\\
A_{i_2}\\
\vdots\\
A_{i_{m'}}
\end{bmatrix}
$$

那么接收端就可以通过解该线性方程组得到 $\vec{c}$ ，**当且仅当 $A'$ 是列满秩的**：

$$
A'\vec{c}=\vec{s'}
$$


#### LT码的编码过程

在本项目中，LT码的编码过程如下：

- 使用伽罗瓦型LFSR（Linear Feedback Shift Register）迭代生成用于为码元采样数据块的随机种子，种子长度为4byte。
- 对于每一个码元，将LFSR生成的随机数作为随机种子，遵循稳定孤波分布生成该码元的度数 $d$ ，并从原始数据中随机采样 $d$ 个数据块来构造负载。在本项目中，数据块的长度为16byte。
- 将LFSR生成的随机种子连接到所有采样到的数据块的异或和之前，构成一个长度为20byte的“水滴”。

不难发现，虽然没有直接记录码元含有哪些数据块，但通过记录随机种子，接收端仍然可以通过使用完全相同的随机分布和随机数生成器来将度数和采样的数据块还原出来。换而言之，可以认为这个随机种子间接地记录了这个码元的“元数据”。

#### LT码的解码算法及其实现

在本项目中，我们采用了后向传播算法，其过程如下：

- 收到一个“水滴”之后，将种子和负载切分开，**根据种子重新生成这个码元的度数和对应的采样到的数据块**。
- 如果某些对应的数据块已知，那么直接将负载与这些已知的数据块进行异或运算，并减去相应的度数。之后：
  - 如果剩余的度数恰好为0，那么负载包含的所有数据块都是已知的，这个“水滴”是冗余的，直接丢弃即可。此时，经过异或运算之后的负载应当为全0。
  - 如果剩余的度数大于1，那么负载包含多个未知的数据块。暂时保留这个负载，记录有哪些数据块和这个负载是对应的。
  - 如果剩余的度数恰好为1，那么负载恰好包含一个目前未知的数据块 $c_{unknown}$ 。直接将这个数据块的值设为负载的值，之后，**对于所有包含这个数据块的、被保留的负载，将这些负载都与 $c_{unknown}$ 进行异或操作，并将度数减一**。通过这样的异或操作，我们就可以将新得到的数据块 $c_{unknown}$ 从这些负载中排除。如果一些负载排除了 $c_{unknown}$ 之后度数变成了1，那么**再在这些负载上进行相同的操作**。这样，我们成功地**将已知数据块的影响传播给了尽可能多的负载**。
- 不断接收新的“水滴”，直到所有的数据块都被恢复。

在`main.py`中，`add_payload`和`propagate`这两个函数是解码算法的核心，`add_payload`接收一个负载并用已知的数据块进行处理，`propagate`尝试将度数恰好为1的负载进行传播。

在`add_payload`中：

```python
    for chunk in chunks:
        if file_chunks[chunk] is not None:
            file_chunk = file_chunks[chunk]
            payload = bytes_xor(payload, file_chunk)
            seen_chunks.add(chunk)
    for chunk in seen_chunks:
        chunks.remove(chunk)
```

这段代码查找已知的数据块，并将负载与已知的数据块进行异或运算。`file_chunks`是一个记录已知数据块的列表，其某一位置如果为`None`则说明该数据块未知，否则为已知。

```python
    for chunk in chunks:
        if chunks_2_payload[chunk] is None:
            chunks_2_payload[chunk] = set()
        chunks_2_payload[chunk].add(payload_idx)
    payload_2_chunks.append(chunks)
    payloads.append(payload)
```

`chunks_2_payload`和`payload_2_chunks`都是元素为`set`的列表，其分别记录了数据块指向含有其的负载的边，和负载指向其含有的数据块的边。这样，我们以**二分图**的形式记录了负载和数据块的对应关系，为传播操作提供了便利。这段代码将新接收的负载的信息记录到二分图，最后，`add_payload`便尝试调用`propagate`来进行传播操作。

在`propagate`中：

```python
    if len(payload_2_chunks[payload_idx]) != 1:
        return
```

如果尝试进行传播的负载不满足“度数为1”，则直接跳过。

```python
    corresponding_chunk_idx = payload_2_chunks[payload_idx].pop()
    file_chunks[corresponding_chunk_idx] = payloads[payload_idx]
    payloads[payload_idx] = bytes_xor(
        payloads[payload_idx],
        file_chunks[corresponding_chunk_idx]
    )  # This should make the payload zero
```

这段代码找出度数为1的负载对应的数据块，并将负载的内容填入原本未知的数据块。

```python
    associated_payloads = chunks_2_payload[corresponding_chunk_idx]
    for associated_payload in associated_payloads:
        # Skip the payload that is propagating
        if associated_payload == payload_idx:
            continue
        payload_2_chunks[associated_payload].remove(corresponding_chunk_idx)
        payloads[associated_payload] = bytes_xor(
            payloads[associated_payload],
            file_chunks[corresponding_chunk_idx]
        )
```

这段代码将所有包含新已知的数据块的负载与新已知的数据块进行异或操作，以此将新已知的数据块从它们中剔除。

```python
    for associated_payload in associated_payloads:
        propagate(associated_payload)
```

最后，对于所有受影响的负载，尝试进行传播。

#### 稳定孤波分布

LT码的度数分布应当被精细设计，因为随意的度数分布可能会带来严重的实用性问题。比如，如果直接将所有码元的度数都设为1，那么解码过程就退化为投球问题（独立、随机地选择集合中的元素，所有元素都被至少选中一次）；对于有 $n$ 个数据块的数据，接收者期望需要接收大约 $n(\ln n+\gamma)$ 个“水滴”才能进行解码，而这是完全无法接受的。另一方面，如果令每个码元都含有较多的数据块，解码工作又会变得困难：不仅异或运算的次数更多，复杂度更高，出现开启传播的度数为1的码元的可能性也会更低，这就违背了“不需要依赖特定码元”的思想。具体地，我们希望能够达到以下易用性要求：**接收端应当能够在接收到尽可能少的码元的情况下成功解码**，以及**在能够成功解码的前提下，度数应尽可能小**。为了应对度数分布的问题，LT码使用了稳定孤波分布，其很好地平衡了解码概率和负载的度数。给定“失败概率” $\delta$ ，以及稳定孤波分布本身需要的参数 $c$ ，设 $R=c\ln\frac{n}{\delta}\sqrt{n}$ ，LT码能够保证接收到约 $n+R\ln\frac{n}{\delta}$ 个“水滴”之后，接收端有不低于 $1-\delta$ 的概率成功解码所有数据块。在本项目中， $\delta$ 的取值为 $0.05$ ， $c$ 的取值为 $0.1$ ，那么在接收到大约 $1904$ 个“水滴”之后，我们就有至少 $0.95$ 的概率解码原文件，与“运行项目”部分的结果相符。

### Python2与Python3的`random`软件包兼容性问题

我们尝试对`robust_soliton.py`在Python2和Python3中的行为进行了比较。我们导出了`PRNG.get_src_blocks_wrap`在Python2和Python3中的输出结果，发现其给出的每个码元的度是相同的，仅仅是采样到的源文件片段不同。基于此，我们猜想`random.sample`函数在Python2和Python3中的表现不同，导致出现问题。

进一步地，我们查看了Python2和Python3中`random.sample`函数的实现，发现Python2和Python3中的`random.sample`函数的实现确实不相同，这验证了我们的猜想。

在Python2中，`random.sample`函数的实现如下：

```python
def sample(self, population, k):
        n = len(population)
        if not 0 <= k <= n:
            raise ValueError("sample larger than population")
        random = self.random
        _int = int
        result = [None] * k
        setsize = 21        # size of a small set minus size of an empty list
        if k > 5:
            setsize += 4 ** _ceil(_log(k * 3, 4)) # table size for big sets
        if n <= setsize or hasattr(population, "keys"):
            # An n-length list is smaller than a k-length set, or this is a
            # mapping type so the other algorithm wouldn't work.
            pool = list(population)
            for i in xrange(k):         # invariant:  non-selected at [0,n-i)
                j = _int(random() * (n-i))
                result[i] = pool[j]
                pool[j] = pool[n-i-1]   # move non-selected item into vacancy
        else:
            try:
                selected = set()
                selected_add = selected.add
                for i in xrange(k):
                    j = _int(random() * n)
                    while j in selected:
                        j = _int(random() * n)
                    selected_add(j)
                    result[i] = population[j]
            except (TypeError, KeyError):   # handle (at least) sets
                if isinstance(population, list):
                    raise
                return self.sample(tuple(population), k)
        return result
```

而在Python3中，`random.sample`函数的实现如下：

```python
    def sample(self, population, k, *, counts=None):
        if not isinstance(population, _Sequence):
            if isinstance(population, _Set):
                _warn('Sampling from a set deprecated\n'
                      'since Python 3.9 and will be removed in a subsequent version.',
                      DeprecationWarning, 2)
                population = tuple(population)
            else:
                raise TypeError("Population must be a sequence.  For dicts or sets, use sorted(d).")
        n = len(population)
        if counts is not None:
            cum_counts = list(_accumulate(counts))
            if len(cum_counts) != n:
                raise ValueError('The number of counts does not match the population')
            total = cum_counts.pop()
            if not isinstance(total, int):
                raise TypeError('Counts must be integers')
            if total <= 0:
                raise ValueError('Total of counts must be greater than zero')
            selections = self.sample(range(total), k=k)
            bisect = _bisect
            return [population[bisect(cum_counts, s)] for s in selections]
        randbelow = self._randbelow
        if not 0 <= k <= n:
            raise ValueError("Sample larger than population or is negative")
        result = [None] * k
        setsize = 21        # size of a small set minus size of an empty list
        if k > 5:
            setsize += 4 ** _ceil(_log(k * 3, 4))  # table size for big sets
        if n <= setsize:
            # An n-length list is smaller than a k-length set.
            # Invariant:  non-selected at pool[0 : n-i]
            pool = list(population)
            for i in range(k):
                j = randbelow(n - i)
                result[i] = pool[j]
                pool[j] = pool[n - i - 1]  # move non-selected item into vacancy
        else:
            selected = set()
            selected_add = selected.add
            for i in range(k):
                j = randbelow(n)
                while j in selected:
                    j = randbelow(n)
                selected_add(j)
                result[i] = population[j]
        return result
```

最终，我们在`robust_soliton.py`中增加了`py2_sample`函数，它是一个符合Python2中的`random.sample`函数行为的随机采样函数实现，由此解决了Python2与Python3的`random`软件包兼容性问题。
