# Developer guide for Phylonco - LPhy and LPhyBEAST extension

## Read first

- [LPhy developer guide](https://github.com/LinguaPhylo/linguaPhylo/blob/master/DEV_NOTE.md)
- [LPhyBEAST developer guide](https://github.com/LinguaPhylo/LPhyBeast/blob/master/DEV_NOTE.md)

## For Phylonco developers

### Project structure

Phylonco project contains 3 subprojects:

1. phylonco-beast (BEAST package)
2. phylonco-lphybeast (LPhyBEAST extension, also BEAST package)
3. phylonco-lphy (LPhy extension)

Please note 1 and 2 will release as one BEAST package, 3 will release as a LPhy extension.

### Gradle build

1. How to update dependencies in Intellij, especially SNAPSHOT version.

https://github.com/LinguaPhylo/LPhyBeast/blob/master/DEV_NOTE.md#update-dependencies-in-intellij

2. How to refresh the LPhyBEAST dependency:

<a href="./InstallLPhyBEAST.png"><img src="InstallLPhyBEAST.png" align="right" height="300" ></a>

LPhyBEAST is released as a `.zip` file into the [Maven repository](https://central.sonatype.com/namespace/io.github.linguaphylo),
because of the requirement of BEAST package framework.
Therefore, it is not straight forwards as other dependencies, when you want to upgrade it to a newer version.
You need to run the task `installLPhyBEAST` inside the build of `phylonco-lphybeast` subproject __twice__,
at the first time it downloads the zip file and unzip it, at the second time it loads all jar files into the system.
This process can be done either from command line below, or IntelliJ (screenshot on the right side).

```bash
./gradlew phylonco-lphybeast:installLPhyBEAST
```


