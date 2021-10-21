How to use git

## Change author for all commits

```
git filter-branch -f --env-filter "GIT_AUTHOR_NAME='Newname'; GIT_AUTHOR_EMAIL='new@email'; GIT_COMMITTER_NAME='Newname'; GIT_COMMITTER_EMAIL='new@email';" HEAD
```
From https://stackoverflow.com/a/750191/3326606

```
git filter-branch -f --env-filter "GIT_AUTHOR_NAME='3ki5tj'; GIT_AUTHOR_EMAIL='3ki5tj@gmail.com'; GIT_COMMITTER_NAME='3ki5tj'; GIT_COMMITTER_EMAIL='3ki5tj@gmail.com';" HEAD
```

For future commits
```
git config --global user.name "3ki5tj"
git config --global user.email "3ki5tj@gmail.com"

```


## List current author name
```
git config -l
```
