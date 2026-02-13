# Developer guide for Phylonco 

### Project structure

Phylonco project contains 3 subprojects:

1. phylonco-beast (BEAST package)
2. phylonco-lphybeast (LPhyBEAST extension, also a BEAST package)
3. phylonco-lphy (LPhy extension)

About releases:
* (1) and (2) will are released as one BEAST package. 
* (3) is released separately as a LPhy extension.

### Quickstart guide for developers

#### Requirements:
* Java 17 or higher with Java FX (we recommend Zulu 17 with FX)
* [Git](https://github.com/git-guides/install-git) or [Gitbash](https://git-scm.com/install/windows) for Windows
* [IntelliJ](https://www.jetbrains.com/idea/download/)

#### Setting up tools and your project
1. To autheticate with github, we will be using SSH keys, see [here](https://docs.github.com/en/enterprise-cloud@latest/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)

2. Next, clone the phylonco project using git ssh
```
git clone git@github.com:bioDS/beast-phylonco.git
```

3. Copy the idea projet settings 
```
cp -r ./IntelliJ/.idea/ .
```

4. Open IntelliJ and clear the [cache](https://www.jetbrains.com/help/idea/invalidate-caches.html)

5. Restart IntelliJ, select `Open` and choose the root directory

6. To build the project in IntelliJ, from the menu go to `Build` and `Rebuild project`

#### SDK and dependencies
The [project settings](https://www.jetbrains.com/help/idea/project-settings-and-structure.html) should show the `SDK: JDK 17` (or similar), and `Language level: 17`. 

When new dependencies are added, you may also need to reload your Maven dependencies in IntelliJ, see [here](https://www.jetbrains.com/help/idea/delegate-build-and-run-actions-to-maven.html#maven_reimport)


### Building releases and jars
1. In IntelliJ, go to the Maven tab on the right. 

2. Under `Profiles` select `skipLPhyTests` 

3. Click on `linguaphylo` then `clean`, `install`

4. Click on `lphybeast.root` then `clean`, `install`

5. Click on `phylonco` then `clean`, `install`

The release (packaged) zip files, and jars will be in a target directory.

### Additional Resources
* [LPhy setup guide](https://github.com/LinguaPhylo/linguaPhylo/blob/master/DEV_NOTE.md)
* [LPhy language guide](https://github.com/LinguaPhylo/linguaPhylo/blob/master/DEV_NOTE2.md)
* [Maven modules](https://github.com/LinguaPhylo/linguaPhylo/blob/master/DEV_NOTE1.md)
* [Maven projects](https://github.com/LinguaPhylo/linguaPhylo/blob/master/DEV_NOTE3.md)
